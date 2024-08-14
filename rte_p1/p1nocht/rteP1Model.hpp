#if !defined(rte_p1_model_)
#define rte_p1_model_

namespace rteP1Model
{
void setup0(const int tid, const int gid, 
            const int twrkid, const int gwrkid,
            const dfloat tau, const dfloat bo, const dfloat eps);
void setup0();
void sendParamToNek();
void buildKernel(occa::properties& kernelInfo);
void setup(nrs_t *nrs, std::vector<int> GbidRobin_);
void extraSetup(nrs_t *nrs);
void solve(nrs_t *nrs, double time, int tstep, occa::memory o_T);
void solve(nrs_t *nrs, double time, int tstep);

void updateProperties(const double time, const int tstep);
void updateForT(nrs_t *nrs, occa::memory o_Tnew, occa::memory o_Gnew);
void checkCHTFlux(nrs_t *nrs, const double time, const int tstep,
                  const dfloat conFluid, const dfloat conSolid,
                  std::vector<int> TbidFluidSide, std::vector<int> TbidSolidSide);

occa::memory o_implicitRobinLHS_g(int scalarIdx);
occa::memory o_radiationSource_t();
}

// private members
namespace
{
static mesh_t *meshT;
static mesh_t *meshG;

// input parameters
int P_RTE_TID;    // scalar id for temperature
int P_RTE_GID;    // scalar id for incident radiation (main solution)
int P_RTE_TWRKID; // index for o_usrwrk
int P_RTE_GWRKID;
dfloat P_RTE_TAU; // optical length
dfloat P_RTE_BO;  // Boltzman number
dfloat P_RTE_EPS; // Emissivity

std::vector<int> GbidRobin; // boundary id for G's Robin BC

// work array
static occa::memory o_GRobinLHS;
static occa::memory o_Glambda0;
static occa::memory o_Glambda1;
//static occa::memory o_Glag;   // TODO: extrapolation?
static occa::memory o_Tsource;
static occa::memory o_invJw;  // TODO: issue with strogGrad

occa::kernel TrhsKernel;
occa::kernel GrhsKernel;
occa::kernel compFluxChtKernel;

dfloat areaCHT_f, areaCHT_s;

// flags
static bool buildKernelCalled = false;
static bool setup0Called = false;
static bool setupCalled = false;
static bool setupHmhSolverCalled = false;
static bool setupComplete = false; // we have too many setup...
static bool checkCHTFluxCalled = false;
} // private members

// default parameters
static dfloat coeff[] = {
   0,         // scalar id of T
   1,         // scalar id of G
   0,         // usrwrk slot id of T
   1,         // usrwrk slot id of G
   4.0,       // tau, optical length
   57.0,      // Bo, Boltzmann number
   0.85,      // eps, emissivity
};

occa::memory rteP1Model::o_implicitRobinLHS_g(int scalarIdx)
{
  if (scalarIdx == P_RTE_GID && GbidRobin.size() > 0) {
    return o_GRobinLHS;
  }
  return o_NULL;
}

occa::memory rteP1Model::o_radiationSource_t()
{
//  occa::memory o_T = nrs->cds->o_S + nrs->cds->fieldOffsetScan[P_RTE_TID];
//  occa::memory o_G = nrs->cds->o_S + nrs->cds->fieldOffsetScan[P_RTE_GID];
//  rteP1Model::updateForT(nrs, o_T, o_G);
  return o_Tsource;
}

// Feed global variables from udf via arguments
void rteP1Model::setup0(const int tid, const int gid,
                        const int twrkid, const int gwrkid, 
                        const dfloat tau, const dfloat bo, const dfloat eps)
{
  nekrsCheck(setup0Called,
             platform->comm.mpiComm,
             EXIT_FAILURE,
             "%s\n",
             "RTE setup0 has already been called!");

  P_RTE_TID = tid;
  P_RTE_GID = gid;
  P_RTE_TWRKID = twrkid;
  P_RTE_GWRKID = gwrkid;

  P_RTE_TAU = tau;
  P_RTE_BO = bo;
  P_RTE_EPS = eps;

  if (platform->comm.mpiRank==0) {
    printf("nekRTE CONST: scalarID for temperature    %d\n", P_RTE_TID);
    printf("nekRTE CONST: scalarID for incident       %d\n", P_RTE_GID);
    printf("nekRTE CONST: usrwrk id for temperature   %d\n", P_RTE_TWRKID);
    printf("nekRTE CONST: usrwrk id for incident      %d\n", P_RTE_GWRKID);
    printf("nekRTE CONST: Optical length              %8.6e\n", P_RTE_TAU);
    printf("nekRTE CONST: Boltzmann number            %8.6e\n", P_RTE_BO);
    printf("nekRTE CONST: Emissivity                  %8.6e\n", P_RTE_EPS);
  }
  fflush(stdout);

  setup0Called = true;
}

// Feed global variables from default or par
void rteP1Model::setup0()
{
  std::string c = "casedata", s;
  int tid    = (platform->par->extract(c, "p1_tid",    s)) ? std::stoi(s) : coeff[0];
  int gid    = (platform->par->extract(c, "p1_gid",    s)) ? std::stoi(s) : coeff[1];
  int twrkid = (platform->par->extract(c, "p1_twrkid", s)) ? std::stoi(s) : coeff[2];
  int gwrkid = (platform->par->extract(c, "p1_gwrkid", s)) ? std::stoi(s) : coeff[3];

  dfloat tau = platform->par->extract(c, "p1_tau", s) ? std::stod(s) : coeff[4];
  dfloat bo  = platform->par->extract(c, "p1_bo",  s) ? std::stod(s) : coeff[5];
  dfloat eps = platform->par->extract(c, "p1_eps", s) ? std::stod(s) : coeff[6];

  rteP1Model::setup0(tid, gid, twrkid, gwrkid, tau, bo, eps);
}

// call this in UDF_Setup0
void rteP1Model::sendParamToNek()
{
  if (platform->options.compareArgs("BUILD ONLY", "FALSE")) {
    *nek::ptr<double>("p_p1_tau") = P_RTE_TAU;
    *nek::ptr<double>("p_p1_bo") = P_RTE_BO;
    *nek::ptr<double>("p_p1_eps") = P_RTE_EPS;
  }
}


void rteP1Model::buildKernel(occa::properties& kernelInfo)
{ 
  nekrsCheck(!setup0Called,
             platform->comm.mpiComm,
             EXIT_FAILURE,
             "%s\n",
             "called prior rteP1Model::setup0()!");

  nekrsCheck(buildKernelCalled,
             platform->comm.mpiComm,
             EXIT_FAILURE,
             "%s\n",
             "RTE buildKernel has already been called!");

  kernelInfo["defines/P_RTE_TID"] = P_RTE_TID;
  kernelInfo["defines/P_RTE_GID"] = P_RTE_GID;
  kernelInfo["defines/P_RTE_TWRKID"] = P_RTE_TWRKID;
  kernelInfo["defines/P_RTE_GWRKID"] = P_RTE_GWRKID;
  kernelInfo["defines/P_RTE_TAU"] = P_RTE_TAU;
  kernelInfo["defines/P_RTE_BO"] = P_RTE_BO;
  kernelInfo["defines/P_RTE_EPS"] = P_RTE_EPS;

  auto buildKernel = [&kernelInfo](const std::string &kernelName) {
    // const auto path = getenv("NEKRS_KERNEL_DIR") + std::string("/nrs/plugins/"); TODO
    const auto path = std::string(fs::current_path()) + "/";
    const auto fileName = path + "rteP1Model.okl";
    const auto reqName = "rteP1Model::";
    if (platform->options.compareArgs("REGISTER ONLY", "TRUE")) {
      platform->kernelRequests.add(reqName, fileName, kernelInfo);
      return occa::kernel();
    } else {
      buildKernelCalled = true;
      return platform->kernelRequests.load(reqName, kernelName);
    }
  };

  TrhsKernel = buildKernel("T_rhs");
  GrhsKernel = buildKernel("G_rhs");
  compFluxChtKernel = buildKernel("compFluxCht");
}

void rteP1Model::updateProperties(const double time, const int tstep)
{
  if (platform->comm.mpiRank == 0) {
    printf("RTE updateProperties %d %g sid=%d",tstep,time,P_RTE_GID);
  }

  { // Robin LHS
    const dfloat hc = P_RTE_EPS / ( 2.0 * P_RTE_TAU * (2.0 - P_RTE_EPS) );

    if (platform->comm.mpiRank == 0) printf(", hc=%12.4e", hc);

    auto o_surfaceAreaRate = [&]() // TODO, only do this once?
    {
      auto nbid = GbidRobin.size();
      auto o_bid = platform->o_memPool.reserve<int>(nbid);
      o_bid.copyFrom(GbidRobin.data(), nbid);

      // compute area / bm1
      auto o_surfaceAreaRate = meshG->surfaceAreaMultiply(nbid, o_bid, o_invJw);
      return o_surfaceAreaRate;
    }();

    platform->linAlg->fill(meshG->Nlocal, hc, o_GRobinLHS);
    platform->linAlg->axmy(meshG->Nlocal, 1.0, o_surfaceAreaRate, o_GRobinLHS); // y[n] = a*x[n]*y[n]
  }

  { // Helmholtz coef
    const dfloat p1_diff = 1.0 / (3.0 * P_RTE_TAU*P_RTE_TAU);
    const dfloat p1_rho = 1.0;

    if (platform->comm.mpiRank == 0) {
      printf(", lambda0=%12.4e, lambda1=%12.4e\n", p1_diff, p1_rho);
    }

    platform->linAlg->fillKernel(meshG->Nlocal, p1_diff, o_Glambda0);
    platform->linAlg->fillKernel(meshG->Nlocal, p1_rho, o_Glambda1);
  }
}

void rteP1Model::updateForT(nrs_t *nrs, occa::memory o_Tnew, occa::memory o_Gnew)
{
  nekrsCheck(!setupCalled,
             platform->comm.mpiComm,
             EXIT_FAILURE,
             "%s\n",
             "called prior rteP1Model::setup()!");

  // fill usrwrk for BC
  nrs->o_usrwrk.copyFrom(o_Tnew, meshT->Nlocal, P_RTE_TWRKID*nrs->fieldOffset);
  nrs->o_usrwrk.copyFrom(o_Gnew, meshG->Nlocal, P_RTE_GWRKID*nrs->fieldOffset);

  // TODO ext?

  // Update Tsrc, T can be cht, but source from G is v-mesh only
  platform->linAlg->fill(meshT->Nlocal, 0.0, o_Tsource);
  TrhsKernel(meshG->Nlocal, o_Tnew, o_Gnew, o_Tsource);
}

void rteP1Model::setup(nrs_t *nrs, std::vector<int> GbidRobin_)
{ 
  auto cds = nrs->cds;
  meshT = (P_RTE_TID) ? nrs->meshV : nrs->cds->mesh[0];
  meshG = nrs->meshV; // G is always on v-mesh

  const int NSfields = (cds) ? cds->NSfields : 0;

  nekrsCheck(!setup0Called,
             platform->comm.mpiComm,
             EXIT_FAILURE,
             "%s\n", 
             "called prior rteP1Model::setup0()!");

  nekrsCheck(!buildKernelCalled,
             platform->comm.mpiComm,
             EXIT_FAILURE,
             "%s\n", 
             "called prior rteP1Model::buildKernel()!");

  nekrsCheck(setupCalled,
             platform->comm.mpiComm,
             EXIT_FAILURE,
             "%s\n", 
             "invalid second call");

  nekrsCheck(NSfields==0,
             platform->comm.mpiComm,
             EXIT_FAILURE,
             "Nscalars cannot be 0, %d\n", 
             NSfields);

  nekrsCheck(P_RTE_TID==NSfields,
             platform->comm.mpiComm,
             EXIT_FAILURE,
             "P_RTE_TID cannot be the same as P_RTE_GID! %d %d\n",
             P_RTE_TID, P_RTE_GID);

  nekrsCheck(P_RTE_TID>=NSfields,
             platform->comm.mpiComm,
             EXIT_FAILURE,
             "P_RTE_TID exceeds nrs->Nscalar, %d\n",
             P_RTE_TID);

  nekrsCheck(P_RTE_GID>=NSfields,
             platform->comm.mpiComm,
             EXIT_FAILURE,
             "P_RTE_GID exceeds nrs->Nscalar, %d\n",
             P_RTE_GID);

  GbidRobin = GbidRobin_;

  if (platform->comm.mpiRank==0) { 
    printf("nekRTE PARAM: Robin bid for G          ");
    for (int id : GbidRobin) printf("%d ", id);
    if (GbidRobin.size()==0) printf("(none)");
    printf("\n");
  }

  // Allocate work arrays
  o_GRobinLHS = platform->device.malloc<dfloat>(cds->fieldOffset[P_RTE_GID]);
  o_Glambda0 = platform->device.malloc<dfloat>(cds->fieldOffset[P_RTE_GID]);
  o_Glambda1 = platform->device.malloc<dfloat>(cds->fieldOffset[P_RTE_GID]);
  o_Tsource = platform->device.malloc<dfloat>(cds->fieldOffset[P_RTE_TID]);

  // TODO: for strongGrad issue
  o_invJw = platform->device.malloc<dfloat>(cds->fieldOffset[P_RTE_TID]);
  o_invJw.copyFrom(meshT->o_LMM, meshT->Nlocal);
  platform->linAlg->ady(meshT->Nlocal, 1.0, o_invJw);

  // Fill variables
  rteP1Model::updateProperties(0.0, 0);

  setupCalled = true;
}

// call at tstep==0
void rteP1Model::extraSetup(nrs_t *nrs) {
  nekrsCheck(!setupCalled,
             platform->comm.mpiComm,
             EXIT_FAILURE,
             "%s\n",  
             "called prior rteP1Model::extraSetup()!");
  if (!setupHmhSolverCalled) {
    hmhSolver::setup(nrs, P_RTE_GID);
    nrs->cds->compute[P_RTE_GID] = 0;
    setupHmhSolverCalled = true;
    setupComplete = true;
  }

  occa::memory o_T = nrs->cds->o_S + nrs->cds->fieldOffsetScan[P_RTE_TID];
  occa::memory o_G = nrs->cds->o_S + nrs->cds->fieldOffsetScan[P_RTE_GID];
  rteP1Model::updateForT(nrs, o_T, o_G);
}

// Interface for user provided Tin, so we can feed Tavg
void rteP1Model::solve(nrs_t *nrs, double time, int tstep, occa::memory o_T)
{
  nekrsCheck(!setupComplete,
             platform->comm.mpiComm,
             EXIT_FAILURE,
             "%s %d %d %d %d\n",  
             "called but some setup is not done!",
             setup0Called, buildKernelCalled, setupCalled, setupHmhSolverCalled);

  occa::memory o_G = nrs->cds->o_S + nrs->cds->fieldOffsetScan[P_RTE_GID];
  occa::memory o_rhs = platform->o_memPool.reserve<dfloat>(meshG->Nlocal);

  // copy for BC
  nrs->o_usrwrk.copyFrom(o_T, meshT->Nlocal, P_RTE_TWRKID*nrs->fieldOffset);
  nrs->o_usrwrk.copyFrom(o_G, meshG->Nlocal, P_RTE_GWRKID*nrs->fieldOffset); // not-used

  GrhsKernel(meshG->Nlocal, o_T, o_rhs); // TODO: extrapolate in time?

  hmhSolver::solve(nrs, P_RTE_GID, time, tstep, o_Glambda0, o_Glambda1, o_rhs);

  rteP1Model::updateForT(nrs, o_T, o_G);
}

// Default, o_T from nrs->o_T
void rteP1Model::solve(nrs_t *nrs, double time, int tstep) 
{
  occa::memory o_T = nrs->cds->o_S + nrs->cds->fieldOffsetScan[P_RTE_TID];
  rteP1Model::solve(nrs, time, tstep, o_T);
}

void rteP1Model::checkCHTFlux(nrs_t *nrs, const double time, const int tstep,
                              const dfloat conFluid, const dfloat conSolid,
                              std::vector<int> TbidFluidSide,
                              std::vector<int> TbidSolidSide)
{
  auto cds = nrs->cds;
  const dlong TfieldOffset = cds->fieldOffset[P_RTE_TID];
  occa::memory o_T = cds->o_S + cds->fieldOffsetScan[P_RTE_TID];
  occa::memory o_G = cds->o_S + cds->fieldOffsetScan[P_RTE_GID];

  // merge two sides
  std::vector<int> bidBoth;
  bidBoth.insert(bidBoth.end(), TbidFluidSide.begin(), TbidFluidSide.end());
  bidBoth.insert(bidBoth.end(), TbidSolidSide.begin(), TbidSolidSide.end());

  // send bid to device
  auto o_copy_bid = [](std::vector<int> bid) 
  {        
    auto nbid = bid.size();
    auto o_bid = platform->o_memPool.reserve<int>(nbid);
    o_bid.copyFrom(bid.data(), nbid);
    return o_bid;
  };
  auto o_bid_f = o_copy_bid(TbidFluidSide);
  auto o_bid_s = o_copy_bid(TbidSolidSide);
  auto o_bidBoth = o_copy_bid(bidBoth);

  // surface area
  if (!checkCHTFluxCalled) {
    auto o_one = platform->o_memPool.reserve<dfloat>(TfieldOffset);
    platform->linAlg->fill(meshT->Nlocal, 1.0, o_one);
    areaCHT_f = meshT->surfaceAreaMultiplyIntegrate(o_bid_f.length(), o_bid_f, o_one).at(0);
    areaCHT_s = meshT->surfaceAreaMultiplyIntegrate(o_bid_s.length(), o_bid_s, o_one).at(0);
    checkCHTFluxCalled = true;
  }

  // compute gradient, TODO: strongGrad issue
  auto my_grad = [](mesh_t *mesh, dlong offset, occa::memory o_fld)
  {
    auto o_grad = opSEM::strongGrad(mesh, offset, o_fld); // weak grad
    platform->linAlg->axmyMany(mesh->Nlocal, 3, offset, 0, 1.0, o_invJw, o_grad);
    return o_grad;
  };
  auto o_gradT = my_grad(meshT, TfieldOffset, o_T);

  // nusselt number
  auto dTdn_f = meshT->surfaceAreaNormalMultiplyIntegrate(
                                      TfieldOffset,
                                      o_bid_f.length(),
                                      o_bid_f,
                                      o_gradT).at(0);

  dfloat nu = dTdn_f / areaCHT_f;

  if (platform->comm.mpiRank == 0) {
    printf("%9d %11.4e %11.4e %11.4e %11.4e nusselt(rs)\n",
           tstep, time, dTdn_f, areaCHT_f, nu);
  }

  if (!nrs->cht) return;

  // compute flux: o_Jdotn = ( - con * grad T - coef*grad G ) dot n
  const dfloat rte_coef = 4.0 / (3.0 * P_RTE_TAU * P_RTE_BO);

  auto o_gradG = my_grad(meshG, TfieldOffset, o_G);
  auto o_Jdotn = platform->o_memPool.reserve<dfloat>(TfieldOffset);
  compFluxChtKernel(meshT->Nelements,
                    TfieldOffset,
                    o_bidBoth.length(),
                    o_bidBoth,
                    meshT->o_sgeo,
                    meshT->o_vmapM,
                    meshT->o_EToB,
                    meshT->o_elementInfo,
                    -conFluid,
                    -conSolid,
                    o_gradT,
                    -rte_coef,
                    o_gradG,
                    o_Jdotn);

  auto flux_f = meshT->surfaceAreaMultiplyIntegrate(o_bid_f.length(), o_bid_f, o_Jdotn).at(0);
  auto flux_s = meshT->surfaceAreaMultiplyIntegrate(o_bid_s.length(), o_bid_s, o_Jdotn).at(0);

  // chk1: total flux
  auto err1 = fabs(flux_f + flux_s);

  // chk2: pointwise flux, err = max abs [[J]] = max abs dssum(Js,Jf);
  oogs::startFinish(o_Jdotn, 1, TfieldOffset, ogsDfloat, ogsAdd, meshT->oogs);
  auto err2 = platform->linAlg->amax(meshT->Nlocal, o_Jdotn, platform->comm.mpiComm);

  // chk3: area chk
  auto err3 = fabs(areaCHT_f - areaCHT_s);

  if (platform->comm.mpiRank==0) {
    printf("rte_cht_flux(rs) %8d %11.4e %11.4e %11.4e ", tstep, time, flux_f, flux_s);
    printf("%11.4e %11.4e %11.4e %11.4e %11.4e\n", areaCHT_f, areaCHT_s, err1, err2, err3);
  }
}
#endif // hpp
