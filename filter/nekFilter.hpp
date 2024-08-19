#if !defined(nek_filter_hpp_)
#define nek_filter__hpp_

// API
namespace nekFilter
{
void buildKernel(occa::properties& kernelInfo);
void setup(std::string tag, mesh_t *mesh, int nModes, double strength); // TODO, vector?
void apply(std::string tag, occa::memory o_fld);
}

struct FilterContainer
{
  int nModes = 0;
  double strength = 0.0;

  mesh_t *mesh;
  occa::memory o_filterRT;

  bool setupCalled = false;
};

// private member
namespace
{
  std::map<std::string, FilterContainer> filterMap;
  occa::kernel applyFilterKernel;
  bool buildKernelCalled = false;
} // private


namespace // local functions copied from src/core/lowPassFilter.cpp
{

double filterFactorial(int n)
{ 
  if (n == 0) {
    return 1;
  } else { 
    return n * filterFactorial(n - 1);
  }
}

// jacobi polynomials at [-1,1] for GLL
double filterJacobiP(double a, double alpha, double beta, int N)
{ 
  double ax = a;
  
  auto P = (double *)calloc((N + 1), sizeof(double));
  
  // Zero order 
  double gamma0 = pow(2, (alpha + beta + 1)) / (alpha + beta + 1) * filterFactorial(alpha) *
                  filterFactorial(beta) / filterFactorial(alpha + beta);
  double p0 = 1.0 / sqrt(gamma0);
  
  if (N == 0) {
    free(P);
    return p0;
  }
  P[0] = p0;
  
  // first order
  double gamma1 = (alpha + 1) * (beta + 1) / (alpha + beta + 3) * gamma0;
  double p1 = ((alpha + beta + 2) * ax / 2 + (alpha - beta) / 2) / sqrt(gamma1);
  if (N == 1) {
    free(P);
    return p1;
  }

  P[1] = p1;

  /// Repeat value in recurrence.
  double aold = 2 / (2 + alpha + beta) * sqrt((alpha + 1.) * (beta + 1.) / (alpha + beta + 3.));
  /// Forward recurrence using the symmetry of the recurrence.
  for (int i = 1; i <= N - 1; ++i) {
    double h1 = 2. * i + alpha + beta;
    double anew =
        2. / (h1 + 2.) *
        sqrt((i + 1.) * (i + 1. + alpha + beta) * (i + 1 + alpha) * (i + 1 + beta) / (h1 + 1) / (h1 + 3));
    double bnew = -(alpha * alpha - beta * beta) / h1 / (h1 + 2);
    P[i + 1] = 1. / anew * (-aold * P[i - 1] + (ax - bnew) * P[i]);
    aold = anew;
  }

  double pN = P[N];
  free(P);
  return pN;
}

void filterVandermonde1D(int N, int Np, double *r, double *V)
{
  int sk = 0;
  for (int i = 0; i <= N; i++) {
    for (int n = 0; n < Np; n++) {
      V[n * Np + sk] = filterJacobiP(r[n], 0, 0, i);
    }
    sk++;
  }
}

void legendre_poly(double *L, const double x, const int N)
{
   L[0] = 1.0;
   L[1] = x;
   for (int j = 2; j <= N; j++) {
      const double dj = j;
      L[j] = ( (2*dj-1) * x * L[j-1] - (j-1) * L[j-2] ) / dj;
   }
}


void filterBubbleFunc1D(int N, int Np, double *r, double *V)
{
  auto Lj = (double *)malloc(Np * sizeof(double));

  for (int j = 0; j <= N; j++) {
    const double z = r[j];
    legendre_poly(Lj, z, N);

    V[0 * Np + j] = Lj[0];
    V[1 * Np + j] = Lj[1];
    for (int i = 2; i < Np; i++) {
      V[i * Np + j] = Lj[i] - Lj[i-2];
    }
  }
  free(Lj);
}
} // namespace func


void nekFilter::buildKernel(occa::properties& kernelInfo)
{ 
  nekrsCheck(buildKernelCalled,
             platform->comm.mpiComm,
             EXIT_FAILURE,
             "%s\n",
             "RTE buildKernel has already been called!");
  
  auto buildKernel = [&kernelInfo](const std::string &kernelName) {
    // const auto path = getenv("NEKRS_KERNEL_DIR") + std::string("/nrs/plugins/"); TODO
    const auto path = std::string(fs::current_path()) + "/";
    const auto fileName = path + "nekFilter.okl";
    const auto reqName = "nekFilter::";
    if (platform->options.compareArgs("REGISTER ONLY", "TRUE")) {
      platform->kernelRequests.add(reqName, fileName, kernelInfo);
      return occa::kernel();
    } else {
      buildKernelCalled = true;
      return platform->kernelRequests.load(reqName, kernelName);
    }
  };
  
  applyFilterKernel = buildKernel("filterRTHex3D");
}

// apply wght to diagonal
void filterFunctionRelaxation1D_new(int Nmodes, int Nc, double wght, double *A) {

  // zero matrix
  for (int n = 0; n < Nmodes*Nmodes; n++) {
    A[n] = 0.0;
  }

  // Set all diagonal to 1
  for (int n = 0; n < Nmodes; n++) {
    A[n * Nmodes + n] = 1.0;
  }

  int k0 = Nmodes - Nc;
  for (int k = k0; k < Nmodes; k++) {
    double amp = wght * ((k + 1.0 - k0) * (k + 1.0 - k0)) / (Nc * Nc);
    A[k + Nmodes * k] = 1.0 - amp;
  }
}

occa::memory lowPassFilterSetup_new(std::string tag, mesh_t *mesh, 
                                    const dlong filterNc, const dfloat filterWght) {

  const int verbose = platform->options.compareArgs("VERBOSE", "TRUE");

  nekrsCheck(filterNc < 1,
             platform->comm.mpiComm,
             EXIT_FAILURE,
             "%s: filterNc must be at least 1, but is set to %d\n",
             tag.c_str(), filterNc);

  nekrsCheck(filterWght < 0 || filterWght > 1,
             platform->comm.mpiComm,
             EXIT_FAILURE,
             "%s: filterWght must in [0,1], but is set to %g\n",
             tag.c_str(), filterWght);

  // Construct Filter Function
  int Nmodes = mesh->N + 1; // N+1, 1D GLL points

  auto V = (double *)calloc(Nmodes * Nmodes, sizeof(double));
  auto A = (double *)calloc(Nmodes * Nmodes, sizeof(double));

  // Construct Filter Function
  filterFunctionRelaxation1D_new(Nmodes, filterNc, filterWght, A); // Main difference

  if (platform->comm.mpiRank == 0) {
    printf("filt trn (rs) %5s:", tag.c_str());
    for (int k = 0; k < Nmodes; k++) {
      printf("%7.4f", A[k + Nmodes * k]);
    }
    printf("\n");
  }

  // Construct Vandermonde Matrix
  {
    auto r = (double *)malloc(mesh->Np * sizeof(double));
    for (int i = 0; i < mesh->Np; i++) {
      r[i] = mesh->r[i];
    }
//  filterVandermonde1D(mesh->N, Nmodes, r, V); // Legrendre Polyn.^T
    filterBubbleFunc1D(mesh->N, Nmodes, r, V);  // bubble of Legendre
    free(r);
  }

  // Invert the Vandermonde
  int INFO;
  int N = Nmodes;
  int LWORK = N * N;
  double *WORK = (double *)calloc(LWORK, sizeof(double));
  int *IPIV = (int *)calloc(Nmodes + 1, sizeof(int));
  double *iV = (double *)calloc(Nmodes * Nmodes, sizeof(double));

  for (int n = 0; n < (Nmodes + 1); n++) {
    IPIV[n] = 1;
  }
  for (int n = 0; n < Nmodes * Nmodes; ++n) {
    iV[n] = V[n];
  }

  dgetrf_(&N, &N, (double *)iV, &N, IPIV, &INFO);
  nekrsCheck(INFO, MPI_COMM_SELF, EXIT_FAILURE, "%s\n", "dgetrf failed");

  dgetri_(&N, (double *)iV, &N, IPIV, (double *)WORK, &LWORK, &INFO);
  nekrsCheck(INFO, MPI_COMM_SELF, EXIT_FAILURE, "%s\n", "dgetri failed");

  // V*A*V^-1 in row major
  char TRANSA = 'N';
  char TRANSB = 'N';
  double ALPHA = 1.0, BETA = 0.0;
  int MD = Nmodes;
  int ND = Nmodes;
  int KD = Nmodes;
  int LDA = Nmodes;
  int LDB = Nmodes;

  double *C = (double *)calloc(Nmodes * Nmodes, sizeof(double));

  int LDC = Nmodes;

  dgemm_(&TRANSA, &TRANSB, &MD, &ND, &KD, &ALPHA, A, &LDA, iV, &LDB, &BETA, C, &LDC);

  TRANSA = 'N';
  TRANSB = 'N';

  dgemm_(&TRANSA, &TRANSB, &MD, &ND, &KD, &ALPHA, V, &LDA, C, &LDB, &BETA, A, &LDC);

  auto o_A = platform->device.malloc<dfloat>(Nmodes * Nmodes);
  {
    auto tmp = (dfloat *)calloc(Nmodes * Nmodes, sizeof(dfloat));
    for (int i = 0; i < Nmodes * Nmodes; i++) {
      tmp[i] = A[i];
    }
    o_A.copyFrom(tmp, o_A.length());
    free(tmp);
  }

  if (verbose && platform->comm.mpiRank == 0) {
    for (int j = 0; j < Nmodes; j++) {
      printf("filt mat (rs) %5s:", tag.c_str());
      for (int i = 0; i < Nmodes; i++) {
        printf("%11.4e", A[i + Nmodes * j]);
      }
      printf("\n");
    }
  }

  free(A);
  free(V);
  free(C);
  free(iV);
  free(IPIV);
  free(WORK);

  return o_A;
}

void nekFilter::setup(std::string tag, mesh_t *mesh_in, int nModes, double strength) {

  nekrsCheck(!buildKernelCalled,
             platform->comm.mpiComm,
             EXIT_FAILURE,
             "%s\n",
             "called prior nekFilter::buildKernel()!");

  FilterContainer m = filterMap[tag];

  nekrsCheck(m.setupCalled,
             platform->comm.mpiComm,
             EXIT_FAILURE,
             "filter %s has already been setup!\n",
             tag.c_str());
  
  m.mesh = mesh_in;
  m.nModes = nModes;
  m.strength = strength;

  if (m.nModes > 0 && m.strength > 0) {
    m.o_filterRT = lowPassFilterSetup_new(tag, mesh_in, nModes, strength);
  }
  m.setupCalled = true;

  filterMap[tag] = m;
}

void nekFilter::apply(std::string tag, occa::memory o_fld) {

  FilterContainer m = filterMap[tag];

  nekrsCheck(!m.setupCalled,
             platform->comm.mpiComm,
             EXIT_FAILURE,
             "filter %s is not setup!\n",
             tag.c_str());

  if (m.nModes > 0 && m.strength > 0) {
    applyFilterKernel(m.mesh->Nelements, m.o_filterRT, o_fld);
  }
}

#endif // hpp
