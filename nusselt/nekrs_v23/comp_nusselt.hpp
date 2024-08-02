#if !defined(comp_nusselt_hpp_)
#define comp_nusselt_hpp_

// API
namespace nusselt
{
void buildKernel(occa::properties kernelInfo);
void setup(mesh_t* mesh_in, const dlong fieldOffset_in);
dfloat compute(nrs_t* nrs, std::vector<int> bIDList, occa::memory o_temperature, const dfloat time, const int tstep, bool print);
}

// private members
namespace
{
mesh_t* mesh;
dlong fieldOffset;
static int Nblock;

static occa::kernel getBCFluxKernel;
static occa::kernel sumReductionKernel;

static occa::memory o_grad;
static occa::memory o_gradX, o_gradY, o_gradZ;

static occa::memory o_flux, o_area;
                                
static dfloat *sum1, *sum2;
static occa::memory o_sum1, o_sum2;
}

void nusselt::buildKernel(occa::properties kernelInfo)
{
  const std::string path = getenv("NEKRS_KERNEL_DIR") + std::string("/plugins/");
                                
  std::string fileName, kernelName;
  const std::string extension = ".okl";
  {                             
    kernelName = "getBCFlux";   
    fileName = path + kernelName + extension;
    getBCFluxKernel = platform->device.buildKernel(fileName, kernelInfo, true);
 
    kernelName = "sumReduction"; 
    fileName = path + kernelName + extension;
    sumReductionKernel = platform->device.buildKernel(fileName, kernelInfo, true);
  }
}

void nusselt::setup(mesh_t* mesh_in, const dlong fieldOffset_in)
{
  mesh = mesh_in;
  fieldOffset = fieldOffset_in;

  const dlong NfpTotal = mesh->Nelements * mesh->Nfaces * mesh->Nfp;
 
  // volume
  o_grad = platform->device.malloc(3 * fieldOffset * sizeof(dfloat));
  o_gradX = o_grad + 0 * fieldOffset * sizeof(dfloat);
  o_gradY = o_grad + 1 * fieldOffset * sizeof(dfloat);
  o_gradZ = o_grad + 2 * fieldOffset * sizeof(dfloat);

  // face
  o_flux = platform->device.malloc(NfpTotal * sizeof(dfloat));
  o_area = platform->device.malloc(NfpTotal * sizeof(dfloat));

  // reduction
  Nblock = (NfpTotal + BLOCKSIZE - 1) / BLOCKSIZE;
  sum1 = (dfloat *)calloc(Nblock, sizeof(dfloat));
  sum2 = (dfloat *)calloc(Nblock, sizeof(dfloat));
  o_sum1 = platform->device.malloc(Nblock * sizeof(dfloat));
  o_sum2 = platform->device.malloc(Nblock * sizeof(dfloat));
}

dfloat nusselt::compute(nrs_t* nrs, std::vector<int> bIDList, occa::memory o_temperature, const dfloat time, const int tstep, bool print)
{
  
  const dlong NfpTotal = mesh->Nelements * mesh->Nfaces * mesh->Nfp;

  dfloat sbuf[2] = {0, 0};
  int iter = 0;
  for (const int &bid : bIDList) { // accumulate flux for multiple bcid

    // TODO, strong grad w/o gs, get rid of nrs
    nrs->gradientVolumeKernel(mesh->Nelements,
                              mesh->o_vgeo,
                              mesh->o_D,
                              fieldOffset,
                              o_temperature,
                              o_grad);

    oogs::startFinish(o_grad, 3, fieldOffset, ogsDfloat, ogsAdd, nrs->gsh);

    platform->linAlg->axmy(mesh->Nlocal, 1.0, mesh->o_invLMM, o_gradX); // y[i] = alpha*x[i]*y[i]
    platform->linAlg->axmy(mesh->Nlocal, 1.0, mesh->o_invLMM, o_gradY);
    platform->linAlg->axmy(mesh->Nlocal, 1.0, mesh->o_invLMM, o_gradZ);
  
    getBCFluxKernel(mesh->Nelements,
                    bid,
                    fieldOffset,
                    o_grad,
                    mesh->o_vmapM,
                    mesh->o_EToB,
                    mesh->o_sgeo,
                    o_area,
                    o_flux);

    // reduction 1: device
    sumReductionKernel(NfpTotal, o_area, o_flux, o_sum1, o_sum2);
    o_sum1.copyTo(sum1, Nblock*sizeof(dfloat)); 
    o_sum2.copyTo(sum2, Nblock*sizeof(dfloat)); 
  
    // reduction 2: host
    if (iter==0) {   // only sum area once
      for (int n = 0; n < Nblock; n++) { 
        sbuf[0] += sum1[n];
      }
    }
    for (int n = 0; n < Nblock; n++) {
       sbuf[1] += sum2[n];
    }
    iter++;
  }

  // reduction 3: MPI
  MPI_Allreduce(MPI_IN_PLACE, sbuf, 2, MPI_DFLOAT, MPI_SUM, platform->comm.mpiComm);

  dfloat nu = sbuf[1]/sbuf[0];
  if (platform->comm.mpiRank == 0 && print) {
    printf("%9d %11.4e %11.4e %11.4e %11.4e nusselt(rs)\n", tstep, time, sbuf[1], sbuf[0], nu);
  }
  return nu;
}
#endif // hpp
