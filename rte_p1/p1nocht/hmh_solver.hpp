#if !defined(hmh_solver_hpp_)
#define hmh_solver_hpp_
#include "nrs.hpp"
#include "avm.hpp"
#include "applyDirichlet.hpp"

// API
//
namespace hmhSolver 
{
void setup(nrs_t *nrs, const int is);
void solve(nrs_t *nrs, const int is, const dfloat timeNew, const int tstep,
           occa::memory o_h1, occa::memory o_h2, occa::memory o_rhs);
}

// private members
namespace
{
std::vector<bool> scalarSetupCalled;
static bool initCalled = false;
std::vector<int> cdsComputeBak;

static std::vector<double> timerTotal; // total time spent in hmhSolver
static int numTotalCalled = 0;
}

static void turnOnCdsCompute(cds_t *cds, const int is_in)
{
  cdsComputeBak.resize(cds->NSfields,0);
  for (int is = 0; is < cds->NSfields; is++) {
    cdsComputeBak.at(is) = cds->compute[is];
    cds->compute[is] = 0;
  }
  cds->compute[is_in] = 1;
}

static void recoverCdsCompute(cds_t *cds)
{
  for (int is = 0; is < cds->NSfields; is++) {
    cds->compute[is] = cdsComputeBak.at(is);
  }
}

static void initialization(const int NSfields)
{
  nekrsCheck(initCalled,
             platform->comm.mpiComm,
             EXIT_FAILURE,
             "%s\n",
             "hmhSolver has already been initialized!");

  nekrsCheck(NSfields==0,
             platform->comm.mpiComm,
             EXIT_FAILURE,
             "%s\n",
             "NSfields cannot be 0!");

  scalarSetupCalled.resize(NSfields, false);
  timerTotal.resize(NSfields, 0);
  platform->timer.addUserStat("myHmhSolver::");
  initCalled = true;
}

// setup elliptic solver
void hmhSolver::setup(nrs_t *nrs, const int is)
{
  nrs = dynamic_cast<nrs_t *>(platform->solver);
  cds_t *cds = nrs->cds;
  const int NSfields = (cds) ? cds->NSfields : 0;

  if (!initCalled) {
    initialization(NSfields);
  }

  nekrsCheck(is >= NSfields,
             platform->comm.mpiComm,
             EXIT_FAILURE,
             "%s\n",
             "scalarId exceed NSfields!");

  nekrsCheck(scalarSetupCalled.at(is),
             platform->comm.mpiComm,
             EXIT_FAILURE,
             "scalarId %d has already been setup!\n",
             is);

  nekrsCheck(!cds->compute[is],
             platform->comm.mpiComm,
             EXIT_FAILURE,
             "scalarId %d needed to be turned on in par!!\n",
             is);

  nekrsCheck(cds->cvodeSolve[is],
             platform->comm.mpiComm,
             EXIT_FAILURE,
             "%s\n",
             "hmhSolver does NOT suppoert CVODE");

  scalarSetupCalled.at(is) = true;
}

// select from nrs_t::printInfo
static void printInfo(cds_t *cds, const int is)
{
  std::string sid = scalarDigitStr(is);

  if (platform->comm.mpiRank == 0) {
    elliptic *solver = cds->solver[is];
    std::string name = "S" + sid;
    const auto [prevProjVecs, nProjVecs] = solver->projectionCounters();
    if (nProjVecs > 0) {
      if (prevProjVecs > 0) {
        printf("%-10s: resNorm0 %.2e  resNorm %.2e  ratio = %.3e  %d/%d\n",
               std::string("proj" + name).c_str(),
               solver->initialResidual(),
               solver->initialGuessResidual(),
               solver->initialResidual() / solver->initialGuessResidual(),
               prevProjVecs,
               nProjVecs);
      }
    }
    printf("%-10s: iter %03d  resNorm0 %.2e  resNorm %.2e\n",
           name.c_str(),
           solver->Niter(),
           solver->initialGuessResidual(),
           solver->finalResidual());
  }
}

void hmhSolver::solve(nrs_t *nrs, const int is, const dfloat timeNew, const int tstep,
                      occa::memory o_h1_in, occa::memory o_h2_in, occa::memory o_rhs_in)
{
  nekrsCheck(!scalarSetupCalled.at(is),
             platform->comm.mpiComm,
             EXIT_FAILURE,
             "%s\n",
             "called prior hmhSolver::setup()!");

  cds_t *cds = nrs->cds;

  const int stage = 1;
  std::string sid = scalarDigitStr(is);
  platform->timer.tic("myHmhSolver::S"+sid);
  auto mesh = (is) ? cds->meshV : cds->mesh[0];

  turnOnCdsCompute(cds,is);

  platform->timer.tic("myHmhSolver::S" + sid + "::setup");

  auto o_diff_i = cds->o_diff.slice(cds->fieldOffsetScan[is], mesh->Nlocal);
  auto o_rho_i = cds->o_rho.slice(cds->fieldOffsetScan[is], mesh->Nlocal);
  { // evaluateProperties
    o_diff_i.copyFrom(o_h1_in, mesh->Nlocal);
    o_rho_i.copyFrom(o_h2_in, mesh->Nlocal);

/*  FIXME: https://github.com/stgeke/nekRS/commit/8cd1cc25d6f0966c34fb0e975e1188248208bf98
    // AVM (from evaluateProperties)
    std::string regularizationMethod;
    platform->options.getArgs("SCALAR" + sid + " REGULARIZATION METHOD", regularizationMethod);
    const bool applyAVM = regularizationMethod.find("AVM_HIGHEST_MODAL_DECAY") != std::string::npos;
    if (applyAVM)
      avm::apply(cds, timeNew, is, cds->o_S);
*/
  }

  // rhs
  auto o_rhs = platform->o_memPool.reserve<dfloat>(cds->fieldOffset[is]);
  o_rhs.copyFrom(o_rhs_in, mesh->Nlocal);
  { // cds->makeNLT
    if (platform->options.compareArgs("SCALAR" + sid + " REGULARIZATION METHOD", "HPFRT")) {
      cds->o_NLT.copyFrom(o_rhs, mesh->Nlocal, cds->fieldOffsetScan[is]);

      cds->filterRTKernel(cds->meshV->Nelements,
                          is,
                          1,
                          cds->o_fieldOffsetScan,
                          cds->o_applyFilterRT,
                          cds->o_filterRT,
                          cds->o_filterS,
                          cds->o_rho,
                          cds->o_S,
                          cds->o_NLT);
      double flops = 6 * mesh->Np * mesh->Nq + 4 * mesh->Np;
      flops *= static_cast<double>(mesh->Nelements);
      platform->flopCounter->add("scalarFilterRT", flops);

      o_rhs.copyFrom(cds->o_NLT, mesh->Nlocal, 0, cds->fieldOffsetScan[is]); // src, count, destOffset, srcOffset
    }
  }
  { // makef
    platform->linAlg->axmy(mesh->Nlocal, 1.0, mesh->o_LMM, o_rhs); // y = a*x*y
    if (platform->verbose) {
      const dfloat debugNorm = platform->linAlg->weightedNorm2Many(mesh->Nlocal,
                                                                   1,
                                                                   cds->fieldOffset[0],
                                                                   mesh->ogs->o_invDegree,
                                                                   o_rhs,
                                                                   platform->comm.mpiComm);
      if (platform->comm.mpiRank == 0) {
        printf("BF scalar=%d norm: %.15e\n", is, debugNorm);
      }
    }
  }

  applyDirichletScalars(nrs, timeNew, nrs->cds->o_S, nrs->cds->o_Se);

  { // cds_t::solve
    cds->neumannBCKernel(mesh->Nelements,
                         1,
                         mesh->o_sgeo,
                         mesh->o_vmapM,
                         mesh->o_EToB,
                         is,
                         timeNew,
                         cds->fieldOffset[is],
                         0,
                         cds->EToBOffset,
                         mesh->o_x,
                         mesh->o_y,
                         mesh->o_z,
                         cds->o_Ue,
                         cds->o_S,
                         cds->o_EToB,
                         cds->o_diff,
                         cds->o_rho,
                         *(cds->o_usrwrk),
                         o_rhs);
  }

  const auto o_lambda0 = o_diff_i;
  const auto o_lambda1 = [&]() // Robin via implicit linear term
  { 
    auto o_l = platform->o_memPool.reserve<dfloat>(mesh->Nlocal);
    if (cds->userImplicitLinearTerm) {
      auto o_implicitLT = cds->userImplicitLinearTerm(timeNew, is);
      if(o_implicitLT.isInitialized()) {
        platform->linAlg
          ->axpbyz(mesh->Nlocal, 1.0, o_rho_i, 1.0, o_implicitLT, o_l);
      } else {
        platform->linAlg->axpby(mesh->Nlocal, 1.0, o_rho_i, 0.0, o_l);
      }
    } else {
      platform->linAlg->axpby(mesh->Nlocal, 1.0, o_rho_i, 0.0, o_l);
    }
    return o_l;
  }();

  // initial guess
  auto o_S = [&]() {
    auto o_S0 = platform->o_memPool.reserve<dfloat>(cds->fieldOffset[is]);
    if (platform->options.compareArgs("SCALAR" + sid + " INITIAL GUESS", "EXTRAPOLATION") && stage == 1) {
      o_S0.copyFrom(cds->o_Se, cds->fieldOffset[is], 0, cds->fieldOffsetScan[is]);
    } else { 
      o_S0.copyFrom(cds->o_S, cds->fieldOffset[is], 0, cds->fieldOffsetScan[is]);
    }
    return o_S0;
  }();
  platform->timer.toc("myHmhSolver::S" + sid + "::setup");
  
  // main solve
  std::string precUpdateBackUp;
  platform->options.getArgs("ELLIPTIC PRECO COEFF FIELD", precUpdateBackUp);
  platform->options.setArgs("ELLIPTIC PRECO COEFF FIELD", "TRUE"); //FIXME gmres

  platform->timer.tic("myHmhSolver::S" + sid + "::solve");
  cds->solver[is]->solve(o_lambda0, o_lambda1, o_rhs, o_S);
  platform->timer.toc("myHmhSolver::S" + sid + "::solve");

  platform->options.setArgs("ELLIPTIC PRECO COEFF FIELD", precUpdateBackUp);
  o_S.copyTo(cds->o_S, cds->fieldOffset[is], cds->fieldOffsetScan[is]);

  recoverCdsCompute(cds);
  platform->timer.toc("myHmhSolver::S"+sid);

  timerTotal.at(is) = platform->timer.query("myHmhSolver::S" + sid, "DEVICE:MAX");
  auto timerTotalCalled = std::accumulate(timerTotal.begin(), timerTotal.end(), 0.0);
  numTotalCalled++;
  platform->timer.set("myHmhSolver::Total", timerTotalCalled, numTotalCalled); // in case of many scalar
  auto tPrec = platform->timer.query("scalar" + sid + " preconditioner", "DEVICE:MAX");
  auto nPrec = platform->timer.count("scalar" + sid + " preconditioner");
  platform->timer.set("myHmhSolver::S" + sid + "::solve::preconditioner", tPrec, nPrec);

  printInfo(cds, is);
}

#endif // hpp

