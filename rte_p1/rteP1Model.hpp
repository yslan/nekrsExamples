#if !defined(rte_p1_model_)
#define rte_p1_model_

namespace rteP1Model
{
void buildKernel(occa::properties& kernelInfo, const int tid_, const int sid_, const int wrkid_, const int Gwrkid_,
                 const dfloat tau_, const dfloat bo_, const dfloat eps_);
void setup(nrs_t *nrs);
void solve(nrs_t *nrs, dfloat time, int tstep);
void updateSourceTerms(nrs_t *nrs);
}

// private members
namespace
{
// kernel constant 
int P_RTE_TID;    // scalar id for temperature
int P_RTE_GID;    // scalar id for incident radiation (main solution)
int P_RTE_TWRKID;   // index for o_usrwrk
int P_RTE_GWRKID;
dfloat P_RTE_TAU; // optical length
dfloat P_RTE_BO;  // Boltzman number
dfloat P_RTE_EPS; // Emissivity

// flags
static bool setupCalled = false;
static bool buildKernelCalled = false;
} // private members

void rteP1Model::buildKernel(occa::properties& kernelInfo,
                             const int tid_, const int sid_, const int Twrkid_, const int Gwrkid_,
                             const dfloat tau_, const dfloat bo_, const dfloat eps_)
{ 

  nekrsCheck(buildKernelCalled,
             platform->comm.mpiComm,
             EXIT_FAILURE,
             "%s\n",
             "RTE buildKernel has already been called!");
  
  P_RTE_TID = tid_;
  P_RTE_GID = sid_;
  P_RTE_TWRKID = Twrkid_;
  P_RTE_GWRKID = Gwrkid_;

  P_RTE_TAU = tau_;
  P_RTE_BO = bo_;
  P_RTE_EPS = eps_;
  
  kernelInfo["defines/P_RTE_TID"] = P_RTE_TID;
  kernelInfo["defines/P_RTE_GID"] = P_RTE_GID;
  kernelInfo["defines/P_RTE_TWRKID"] = P_RTE_TWRKID;
  kernelInfo["defines/P_RTE_GWRKID"] = P_RTE_GWRKID;
  kernelInfo["defines/P_RTE_TAU"] = P_RTE_TAU;
  kernelInfo["defines/P_RTE_BO"] = P_RTE_BO;
  kernelInfo["defines/P_RTE_EPS"] = P_RTE_EPS;
  
  if (platform->comm.mpiRank==0) {
    printf("nekRTE CONST: scalarID for temperature    %d\n", P_RTE_TID);
    printf("nekRTE CONST: scalarID for incident       %d\n", P_RTE_GID);
    printf("nekRTE CONST: usrwrk id for temperature   %d\n", P_RTE_TWRKID);
    printf("nekRTE CONST: usrwrk id for incident      %d\n", P_RTE_GWRKID);
    printf("nekRTE CONST: optical length              %8.6e\n", P_RTE_TAU);
    printf("nekRTE CONST: Boltzmann number            %8.6e\n", P_RTE_BO);
    printf("nekRTE CONST: Emissivity                  %8.6e\n", P_RTE_EPS);
  }
  fflush(stdout);

  // UDF_LoadKernels will be called twice
  if (!platform->options.compareArgs("REGISTER ONLY", "TRUE")) {
    buildKernelCalled = true;
  }
}

void rteP1Model::setup(nrs_t *nrs)
{ 
  cds_t *cds = nrs->cds;
  const int NSfields = (cds) ? cds->NSfields : 0;

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

  //nrs = dynamic_cast<nrs_t *>(platform->solver);
  hmhSolver::setup(nrs, P_RTE_GID); 

  setupCalled = true;
}

void rteP1Model::solve(nrs_t *nrs, double time, int tstep) {
  nekrsCheck(!setupCalled,
             platform->comm.mpiComm,
             EXIT_FAILURE,
             "%s\n",  
             "called prior rteP1Model::setup()!");

  cds_t *cds = nrs->cds;
  mesh_t *mesh = nrs->meshV;
  occa::memory o_T = cds->o_S + cds->fieldOffsetScan[P_RTE_TID];
  occa::memory o_G = cds->o_S + cds->fieldOffsetScan[P_RTE_GID];
  
  occa::memory o_h1 = platform->o_memPool.reserve<dfloat>(mesh->Nlocal);
  occa::memory o_h2 = platform->o_memPool.reserve<dfloat>(mesh->Nlocal);
  occa::memory o_rhs = platform->o_memPool.reserve<dfloat>(mesh->Nlocal);
  
  const dfloat p1_diff = 1.0 / (3.0 * P_RTE_TAU*P_RTE_TAU);
  const dfloat p1_rho = 1.0;

  platform->linAlg->fillKernel(mesh->Nlocal, p1_diff, o_h1);
  platform->linAlg->fillKernel(mesh->Nlocal, p1_rho, o_h2);
  
  rteP1Model_compute_P1_rhs(mesh->Nlocal, o_T, o_rhs); // TODO: extrapolate in time?

  hmhSolver::solve(nrs, P_RTE_GID, time, tstep, o_h1, o_h2, o_rhs);
}

void rteP1Model::updateSourceTerms(nrs_t *nrs)
{
  nekrsCheck(!setupCalled,
             platform->comm.mpiComm,
             EXIT_FAILURE,
             "%s\n",
             "called prior rteP1Model::setup()!");
  // TODO ext?

  cds_t *cds = nrs->cds;  
  mesh_t *mesh = nrs->meshV; // T can be cht, but source is fluid-only
  occa::memory o_T = cds->o_S + cds->fieldOffsetScan[P_RTE_TID];
  occa::memory o_G = cds->o_S + cds->fieldOffsetScan[P_RTE_GID];
  occa::memory o_FS = cds->o_NLT + cds->fieldOffsetScan[P_RTE_TID];

  rteP1Model_compute_T_rhs(mesh->Nlocal, o_T, o_G, o_FS);
}
#endif // hpp
