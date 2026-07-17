// fieldExtract.hpp — NekRS v26 box data sampler (header-only)
//
// Interpolates solution fields onto a uniform 1D/2D/3D grid of points (via
// pointInterpolation_t / findpts) and writes one VTK .vts StructuredGrid file
// per call through MPI-IO.
//
// Usage (in a .udf): construct once in UDF_Setup, then call process(time,tstep)
// on each dump step. If the sampler is a global unique_ptr, release it on the
// last step (before MPI_Finalize) so the findpts MPI comm is freed while MPI is
// still live:  if (nrs->lastStep) { sampler.reset(); }
//
//   // Box mode: uniform grid spanning the diagonal x0 -> x1
//   fieldExtract(mesh, boxDims, fldList, fname, x0, x1);
//
//   // Points mode: user-supplied coordinates (one vector per axis)
//   fieldExtract(mesh, boxDims, fldList, fname, XYZ);
//
//   boxDims : point counts {nx, ny, nz}, endpoints included (1D/2D via ny/nz=1)
//   fldList : vector of {name, {deviceMemory<dfloat> per component}}
//   fname   : output prefix -> <fname>{Line|Plane|Box}NNNNN.vts
//
// Average mode (box only): reduce the current fields over one or more axes and
// dump one .vts. Implemented as a custom integer-id gather-scatter (gslib/oogs,
// same pattern as NekRS planarAvg): local points pre-reduce into one partial
// sum per touched reduced cell, ogsAdd combines the cells across ranks (ids =
// reduced-cell index), and the result scatters back. One handle is cached per
// direction and freed in the destructor.
//   sampler->doAvg(time, tstep, "xy");                  // gather: collapse x,y -> line
//   sampler->doAvg(time, tstep, "z", "gather-scatter"); // broadcast z-avg onto the box
//   -> <fname>_avg<dir>{Point|Line|Plane|Box}NNNNN.vts
// It is an integral (trapezoidal) average over the uniform SAMPLING grid, not a
// mass-matrix average over the SEM mesh.
//
// Grid points follow lexicographic (nz, ny, nx) ordering and are distributed
// evenly across MPI ranks. The .vts layout is consumed by test.py / fe_read.py.
//
// Roadmap (see plans/): cylinder / Gauss-Radau layouts, derived fields (dT/dz),
// redistribution of mismatched explicit points, in-memory data() accessor.

#if !defined(nekrs_field_extract_hpp_)
#define nekrs_field_extract_hpp_

#include <map>
#include <memory>

#include "platform.hpp"
#include "mesh3D.h"
#include "pointInterpolation.hpp"

class fieldExtract
{
public:
  using field = std::tuple<std::string, std::vector<deviceMemory<dfloat>>>;

  // Box mode: uniform grid spanning the diagonal x0 -> x1.
  fieldExtract(mesh_t *mesh,
               const std::vector<int> &boxDims,
               const std::vector<field> &fldList,
               const std::string fileName,
               const std::array<dfloat, 3> &x0,
               const std::array<dfloat, 3> &x1);

  // Points mode: user-supplied coordinates (one vector per axis).
  fieldExtract(mesh_t *mesh,
               const std::vector<int> &boxDims,
               const std::vector<field> &fldList,
               const std::string fileName,
               const std::array<std::vector<dfloat>, 3> &XYZ);

  ~fieldExtract(); // frees the cached gather-scatter handles (avgGsh)

  void process(const double time, const int tstep); // interpolate + write .vts

  // Average current fields over avgDir ("x".."xyz", order-insensitive) and dump
  // one .vts. mode: "gather" (collapse averaged axes) | "gather-scatter"
  // (broadcast the average back onto the box). Box mode only.
  void
  doAvg(const double time, const int tstep, const std::string &avgDir, const std::string &mode = "gather");

  // Even point split used internally; exposed so callers building explicit XYZ
  // (points mode) can produce exactly this rank's local share of numPoints.
  static void pointDistribution(dlong numPoints, dlong &numLocal, dlong &offset);

private:
  std::unique_ptr<pointInterpolation_t> interpolator;

  std::vector<field> userFieldList;
  int Nfields;

  std::string fileName;
  int stepCounter;

  int boxDim = 0;
  std::vector<int> boxDims;

  std::array<dfloat, 3> point0{0, 0, 0}; // start pt (box mode; .vts tag only)
  std::array<dfloat, 3> point1{0, 0, 0}; // end pt

  dlong numPoints;      // global number of sample points
  dlong numPointsLocal; // points owned by this rank
  dlong numPointsScan;  // global offset of this rank's first point

  std::vector<dfloat> xCoord, yCoord, zCoord;
  std::vector<size_t> fldDataOffsetScan; // per-field offset into fldData

  std::vector<float> fldData;

  bool boxMode = false; // true: box ctor (x0/x1 grid); false: points ctor
  bool setupCalled = false;
  bool setPointsCalled = false;

  std::map<std::string, int> avgCounter; // per-stream dump counter for doAvg

  // Cached gather-scatter machinery for one doAvg configuration: local box
  // points pre-reduce into per-cell "slot" partial sums (pointSlot maps each
  // local point to its slot), the oogs exchange runs over the slots plus this
  // rank's share of the reduced grid, and both modes read from the one result
  // (gather-scatter: back through pointSlot; gather: the reduced-grid tail).
  struct avgGs_t {
    oogs_t *gsh = nullptr;
    dlong nSlots = 0;             // locally-touched reduced cells
    dlong nTotal = 0;             // nSlots + rLocal (gs vector length)
    dlong rLocal = 0;             // this rank's share of the reduced grid
    dlong rScan = 0;              // global offset of that share
    std::vector<dlong> pointSlot; // local box point -> slot
  };

  std::map<std::string, avgGs_t> avgGsh; // one entry per avg configuration

  void getIxyz(const dlong i, const int nx, const int nxy, int &ix, int &iy, int &iz);

  // Build (or fetch the cached) gslib/oogs machinery for one avg direction.
  avgGs_t &setupAvgGs(const std::array<bool, 3> &avg, const std::string &dirTag);
  void setBoxDims(const std::vector<int> &boxDims_in);

  void setupCommon(const std::vector<int> &boxDims_in,
                   const std::vector<field> &fldList,
                   const std::string &fileName_in);
  void finishSetup(mesh_t *mesh);

  void setPoints(const std::array<dfloat, 3> &x0, const std::array<dfloat, 3> &x1);
  void setPoints(const std::array<std::vector<dfloat>, 3> &XYZ);

  void interpolate();

  // A grid + its interpolated data to emit as one .vts (full box or averaged).
  struct GridView {
    std::array<int, 3> dims; // nx, ny, nz (1 for collapsed axes)
    dlong nGlobal;
    dlong nLocal;
    dlong nScan;
    const std::vector<dfloat> &x;
    const std::vector<dfloat> &y;
    const std::vector<dfloat> &z;
    const std::vector<float> &data;
    const std::vector<size_t> &offsetScan;
    std::array<dfloat, 3> p0;
    std::array<dfloat, 3> p1;
  };

  void writeVts(const GridView &v, const std::string &fname, const double time, const int tstep);
};

inline void fieldExtract::getIxyz(const dlong i, const int nx, const int nxy, int &ix, int &iy, int &iz)
{
  iz = i / nxy;
  int rem = i % nxy;
  iy = rem / nx;
  ix = rem % nx;
}

inline void fieldExtract::setBoxDims(const std::vector<int> &boxDims_in)
{
  boxDims = boxDims_in;
  boxDim = boxDims_in.size();

  nekrsCheck(boxDim < 1 || boxDim > 3,
             MPI_COMM_SELF,
             EXIT_FAILURE,
             "invalid dimensions of the box %d! Support 1d/2d/3d\n",
             boxDim);

  numPoints = 1;
  for (int npt : boxDims) {
    numPoints *= npt;
  }
  nekrsCheck(numPoints < 1, MPI_COMM_SELF, EXIT_FAILURE, "invalid number of points: %d!\n", numPoints);

  pointDistribution(numPoints, numPointsLocal, numPointsScan);
}

inline void fieldExtract::pointDistribution(dlong numPoints, dlong &numLocal, dlong &offset)
{
  const int rank = platform->comm.mpiRank();
  const int Nrank = platform->comm.mpiCommSize();
  const dlong base = numPoints / Nrank; // first evenly distribute points
  const dlong left = numPoints % Nrank; // assign remaining points
  numLocal = (rank < left) ? base + 1 : base;
  offset = (rank > 0) ? rank * base + std::min(left, static_cast<dlong>(rank)) : 0;
}

inline void fieldExtract::setPoints(const std::array<dfloat, 3> &x0, const std::array<dfloat, 3> &x1)
{
  // uniform grid from x0 to x1
  if (numPointsLocal) {
    const int nx = boxDims[0];
    const int ny = (boxDim > 1) ? boxDims[1] : 1;
    const int nz = (boxDim > 2) ? boxDims[2] : 1;
    const int nxy = nx * ny;

    const dfloat dx = (nx > 1) ? (x1[0] - x0[0]) / (nx - 1) : 0.0;
    const dfloat dy = (ny > 1) ? (x1[1] - x0[1]) / (ny - 1) : 0.0;
    const dfloat dz = (nz > 1) ? (x1[2] - x0[2]) / (nz - 1) : 0.0;

    xCoord.resize(numPointsLocal);
    yCoord.resize(numPointsLocal);
    zCoord.resize(numPointsLocal);

    for (dlong i = 0; i < numPointsLocal; i++) {
      const dlong ig = i + numPointsScan;
      int ix, iy, iz;
      getIxyz(ig, nx, nxy, ix, iy, iz);
      xCoord[i] = x0[0] + ix * dx;
      yCoord[i] = x0[1] + iy * dy;
      zCoord[i] = x0[2] + iz * dz;
    }
  }

  setPointsCalled = true;
}

inline void fieldExtract::setPoints(const std::array<std::vector<dfloat>, 3> &XYZ)
{
  const dlong Nlocal = XYZ[0].size();

  { // chk XYZ: all three axes must have matching length
    int err = 0;
    for (int d = 1; d < 3; d++) {
      if (XYZ[d].size() != XYZ[0].size()) {
        err = d;
        break;
      }
    }

    MPI_Allreduce(MPI_IN_PLACE, &err, 1, MPI_INT, MPI_MAX, platform->comm.mpiComm()); // for err-print

    nekrsCheck(err != 0, MPI_COMM_SELF, EXIT_FAILURE, "XYZ[%d] has different length as XYZ[0]!\n", err);
  }

  {                                      // chk XYZ vs boxDims
    nekrsCheck(Nlocal != numPointsLocal, // TODO(later): redistribute points
               MPI_COMM_SELF,
               EXIT_FAILURE,
               "input points don't match the boxDims! Nlocal=%d  numPointsLocal=%d\n",
               Nlocal,
               numPointsLocal);

    int Ntotal = Nlocal;
    MPI_Allreduce(MPI_IN_PLACE, &Ntotal, 1, MPI_INT, MPI_SUM, platform->comm.mpiComm());

    nekrsCheck(Ntotal != numPoints,
               MPI_COMM_SELF,
               EXIT_FAILURE,
               "input points don't match the boxDims! Ntotal=%d  numPoints=%d\n",
               Ntotal,
               numPoints);
  }

  xCoord.resize(numPointsLocal);
  yCoord.resize(numPointsLocal);
  zCoord.resize(numPointsLocal);

  for (dlong i = 0; i < numPointsLocal; ++i) {
    xCoord[i] = XYZ[0][i];
    yCoord[i] = XYZ[1][i];
    zCoord[i] = XYZ[2][i];
  }

  setPointsCalled = true;
}

inline void fieldExtract::setupCommon(const std::vector<int> &boxDims_in,
                                      const std::vector<field> &fldList,
                                      const std::string &fileName_in)
{
  static bool statAdded = false;
  if (!statAdded) {
    platform->timer.addUserStat("fieldExtract::");
    statAdded = true;
  }

  nekrsCheck(setupCalled, MPI_COMM_SELF, EXIT_FAILURE, "%s\n", "Attempt to call setup twice!\n");

  setBoxDims(boxDims_in);

  userFieldList = fldList;
  Nfields = userFieldList.size();

  // Each field's components must share a common length.
  int ifld = 1;
  for (auto &entry : userFieldList) {
    const auto &o_fld = std::get<1>(entry);
    const auto dim_fld = o_fld.size();
    const auto len_fld = o_fld[0].size();

    int err = 0;
    for (size_t idim = 1; idim < dim_fld; idim++) {
      if (o_fld[idim].size() != len_fld) {
        err = idim;
        break;
      }
    }

    MPI_Allreduce(MPI_IN_PLACE, &err, 1, MPI_INT, MPI_MAX, platform->comm.mpiComm()); // for err-print
    nekrsCheck(err != 0,
               MPI_COMM_SELF,
               EXIT_FAILURE,
               "userFieldList[%d][%d] has different length than its 0th array!\n",
               ifld,
               err);
    ifld++;
  }

  fileName = fileName_in;
  stepCounter = 0;
}

inline void fieldExtract::finishSetup(mesh_t *mesh)
{
  nekrsCheck(!setPointsCalled, MPI_COMM_SELF, EXIT_FAILURE, "%s\n", "points are not setup properly!");

  // Layout of interpolated data in fldData: one contiguous block per field,
  // each block laid out component-major (d) then point (i).
  fldDataOffsetScan.resize(Nfields + 1);
  fldDataOffsetScan[0] = 0;
  int ifld = 1;
  for (auto &entry : userFieldList) {
    const auto dim_fld = std::get<1>(entry).size();
    fldDataOffsetScan[ifld] = fldDataOffsetScan[ifld - 1] + numPointsLocal * dim_fld;
    ifld++;
  }

  fldData.resize(fldDataOffsetScan[Nfields]);

  interpolator = std::make_unique<pointInterpolation_t>(mesh, platform->comm.mpiComm());
  interpolator->setPoints(xCoord, yCoord, zCoord);
  interpolator->find();

  setupCalled = true;
}

inline fieldExtract::fieldExtract(mesh_t *mesh,
                                  const std::vector<int> &boxDims_in,
                                  const std::vector<field> &fldList,
                                  const std::string fileName_in,
                                  const std::array<dfloat, 3> &x0,
                                  const std::array<dfloat, 3> &x1)
{
  boxMode = true;
  setupCommon(boxDims_in, fldList, fileName_in);
  point0 = x0;
  point1 = x1;
  setPoints(x0, x1);
  finishSetup(mesh);
}

inline fieldExtract::fieldExtract(mesh_t *mesh,
                                  const std::vector<int> &boxDims_in,
                                  const std::vector<field> &fldList,
                                  const std::string fileName_in,
                                  const std::array<std::vector<dfloat>, 3> &XYZ)
{
  boxMode = false;
  setupCommon(boxDims_in, fldList, fileName_in);
  // point0/point1 stay {0,0,0}; they only feed the .vts x0/x1 tags.
  setPoints(XYZ);
  finishSetup(mesh);
}

inline fieldExtract::~fieldExtract()
{
  // Free every cached doAvg handle (one per avg configuration). Do NOT use
  // oogs::destroy(): it releases the array-new'd oogs_t with free(), which
  // glibc rejects ("free(): invalid pointer" -- the new[] array cookie offsets
  // the pointer), and it does not free the underlying ogs_t either. Instead
  // free the ogs_t via ogsFree and delete[] the handle; the occa::memory
  // members release their buffers in their destructors.
  for (auto &entry : avgGsh) {
    if (entry.second.gsh) {
      ogsFree(entry.second.gsh->ogs);
      delete[] entry.second.gsh;
      entry.second.gsh = nullptr;
    }
  }
  avgGsh.clear();
}

inline fieldExtract::avgGs_t &fieldExtract::setupAvgGs(const std::array<bool, 3> &avg,
                                                       const std::string &dirTag)
{
  auto it = avgGsh.find(dirTag);
  if (it != avgGsh.end()) {
    return it->second;
  }

  const int nx = boxDims[0];
  const int ny = (boxDim > 1) ? boxDims[1] : 1;
  const int nz = (boxDim > 2) ? boxDims[2] : 1;
  const int nxy = nx * ny;

  // reduced grid: averaged axes collapse to 1
  const int rn[3] = {avg[0] ? 1 : nx, avg[1] ? 1 : ny, avg[2] ? 1 : nz};
  const dlong rN = static_cast<dlong>(rn[0]) * rn[1] * rn[2];
  const int rnxy = rn[0] * rn[1];

  avgGs_t ags;
  pointDistribution(rN, ags.rLocal, ags.rScan);

  // Local pre-reduction slots: one per reduced cell touched by this rank's box
  // points. The gs then runs over partial sums (<= 2 local entries per cell:
  // the slot and, if owned here, the reduced-grid tail entry), which keeps the
  // local rowsize tiny no matter how many points collapse into a cell -- a
  // direct point-level gs would exceed ogsSetup's per-row device blocksize
  // limit (e.g. avgDir "xyz" puts the whole box into one row).
  std::map<dlong, dlong> cellSlot; // reduced cell -> slot
  ags.pointSlot.resize(numPointsLocal);
  for (dlong i = 0; i < numPointsLocal; i++) {
    const dlong ig = i + numPointsScan;
    int ix, iy, iz;
    getIxyz(ig, nx, nxy, ix, iy, iz);
    const int rx = avg[0] ? 0 : ix, ry = avg[1] ? 0 : iy, rz = avg[2] ? 0 : iz;
    const dlong rcell = static_cast<dlong>(rz) * rnxy + static_cast<dlong>(ry) * rn[0] + rx;
    const auto ins = cellSlot.emplace(rcell, static_cast<dlong>(cellSlot.size()));
    ags.pointSlot[i] = ins.first->second;
  }
  ags.nSlots = cellSlot.size();
  ags.nTotal = ags.nSlots + ags.rLocal;

  // ids (1-based): slots carry their reduced-cell id; the tail entries carry
  // the id of the reduced point this rank owns (they contribute zero to the
  // ogsAdd sum and receive it back -- that is the gather-mode output)
  std::vector<long long int> ids(ags.nTotal);
  for (const auto &cs : cellSlot) {
    ids[cs.second] = 1 + static_cast<long long int>(cs.first);
  }
  for (dlong i = 0; i < ags.rLocal; i++) {
    ids[ags.nSlots + i] = 1 + static_cast<long long int>(ags.rScan) + i;
  }

  int nComp = 0; // total number of scalar components over all fields
  for (auto &entry : userFieldList) {
    nComp += std::get<1>(entry).size();
  }

  auto ogsh = ogsSetup(ags.nTotal, ids.data(), platform->comm.mpiComm(), 0, platform->device.occaDevice());
  ags.gsh =
      oogs::setup(ogsh, nComp, std::max(ags.nTotal, static_cast<dlong>(1)), ogsDfloat, nullptr, OOGS_AUTO);

  return avgGsh.emplace(dirTag, std::move(ags)).first->second;
}

inline void fieldExtract::interpolate()
{
  platform->timer.tic("fieldExtractnek::interpolate");
  auto o_wrk = platform->deviceMemoryPool.reserve<dfloat>(numPointsLocal);
  std::vector<dfloat> wrk(numPointsLocal);

  int ifld = 0;
  for (const auto &entry : userFieldList) {
    const auto &fieldVector = std::get<1>(entry);
    const auto fieldDim = fieldVector.size();

    for (int d = 0; d < fieldDim; d++) {
      interpolator->eval(1, 0, fieldVector[d], numPointsLocal, o_wrk);
      o_wrk.copyTo(wrk.data(), numPointsLocal);

      for (dlong i = 0; i < numPointsLocal; i++) {
        fldData[fldDataOffsetScan[ifld] + d * numPointsLocal + i] = static_cast<float>(wrk[i]);
      }
    }
    ifld++;
  }
  platform->timer.toc("fieldExtractnek::interpolate");
}

inline void
fieldExtract::writeVts(const GridView &v, const std::string &fname, const double time, const int tstep)
{
  MPI_Comm mpi_comm = platform->comm.mpiComm();
  const int rank = platform->comm.mpiRank();

  // Local point offset in the global point ordering.
  const long long int pOffset = v.nScan;

  // Remove old file / truncate cleanly before MPI-IO opens it.
  if (rank == 0) {
    std::ofstream file(fname, std::ios::binary | std::ios::trunc);
    file.close();
  }
  MPI_Barrier(mpi_comm);

  MPI_File file_out;
  MPI_File_open(mpi_comm, fname.c_str(), MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &file_out);

  auto indent = [](const int n) { return std::string(2 * n, ' '); };

  auto fieldDataArrayInt3 = [](const std::string &name, const int a, const int b, const int c) {
    return "<DataArray type=\"Int32\" Name=\"" + name +
           "\" NumberOfComponents=\"3\" NumberOfTuples=\"1\" format=\"ascii\"> " + std::to_string(a) + " " +
           std::to_string(b) + " " + std::to_string(c) + " </DataArray>\n";
  };

  auto fieldDataArrayFloat3 = [](const std::string &name, const double a, const double b, const double c) {
    return "<DataArray type=\"Float64\" Name=\"" + name +
           "\" NumberOfComponents=\"3\" NumberOfTuples=\"1\" format=\"ascii\"> " + std::to_string(a) + " " +
           std::to_string(b) + " " + std::to_string(c) + " </DataArray>\n";
  };

  auto fieldDataArrayInt64 = [](const std::string &name, const long long v) {
    return "<DataArray type=\"Int64\" Name=\"" + name + "\" NumberOfTuples=\"1\" format=\"ascii\"> " +
           std::to_string(v) + " </DataArray>\n";
  };

  auto vtsDataArray = [](const std::string &fieldName, long long int nComponent, long long int offset) {
    return "<DataArray type=\"Float32\" Name=\"" + fieldName + "\" NumberOfComponents=\"" +
           std::to_string(nComponent) + "\" format=\"appended\" offset=\"" + std::to_string(offset) +
           "\"/>\n";
  };

  const int nx = v.dims[0];
  const int ny = v.dims[1];
  const int nz = v.dims[2];
  const std::string vtsExtent = "0 " + std::to_string(nx - 1) + " " + "0 " + std::to_string(ny - 1) + " " +
                                "0 " + std::to_string(nz - 1);

  constexpr long long int dim = 3;
  constexpr long long int headerBytes = sizeof(unsigned long long int);

  // ---------------------------------------------------------------------------
  // Compute VTK appended-data offsets.
  //
  // These offsets are relative to the byte immediately after the "_" marker.
  // Each appended block is:
  //
  //   [UInt64 nbytes][raw float data]
  //
  // ---------------------------------------------------------------------------
  long long int vtsOffset = 0;

  const long long int posVtsOffset = vtsOffset;
  vtsOffset += headerBytes + dim * v.nGlobal * sizeof(float);

  std::vector<long long int> fieldVtsOffsets;
  fieldVtsOffsets.reserve(Nfields);

  for (const auto &entry : userFieldList) {
    const auto &fieldVector = std::get<1>(entry);
    const int fieldDim = fieldVector.size();

    fieldVtsOffsets.push_back(vtsOffset);
    vtsOffset += headerBytes + fieldDim * v.nGlobal * sizeof(float);
  }

  const long long int totalAppendedBytes = vtsOffset;

  // ---------------------------------------------------------------------------
  // Write XML header.
  // ---------------------------------------------------------------------------
  std::string message;
  message += indent(0) + "<VTKFile type=\"StructuredGrid\" version=\"1.0\" "
                         "byte_order=\"LittleEndian\" header_type=\"UInt64\">\n";
  message += indent(1) + "<StructuredGrid WholeExtent=\"" + vtsExtent + "\">\n";

  message += indent(2) + "<FieldData>\n";
  message += indent(3) +
             "<DataArray type=\"Float32\" Name=\"TimeValue\" "
             "NumberOfTuples=\"1\" format=\"ascii\"> " +
             std::to_string(time) + " </DataArray>\n";
  message += indent(3) +
             "<DataArray type=\"Int32\" Name=\"CYCLE\" "
             "NumberOfTuples=\"1\" format=\"ascii\"> " +
             std::to_string(tstep) + " </DataArray>\n";
  message += indent(3) + fieldDataArrayInt3("boxDims", nx, ny, nz);
  message += indent(3) + fieldDataArrayInt64("numPoints", static_cast<long long>(v.nGlobal));
  message += indent(3) + fieldDataArrayFloat3("x0", v.p0[0], v.p0[1], v.p0[2]);
  message += indent(3) + fieldDataArrayFloat3("x1", v.p1[0], v.p1[1], v.p1[2]);
  message += indent(2) + "</FieldData>\n";

  message += indent(2) + "<Piece Extent=\"" + vtsExtent + "\">\n";
  message += indent(3) + "<Points>\n";
  message += indent(4) + vtsDataArray("Position", dim, posVtsOffset);
  message += indent(3) + "</Points>\n";
  message += indent(3) + "<PointData>\n";

  int ifld = 0;
  for (const auto &entry : userFieldList) {
    const std::string fieldName = std::get<0>(entry);
    const auto &fieldVector = std::get<1>(entry);
    const int fieldDim = fieldVector.size();
    message += indent(4) + vtsDataArray(fieldName, fieldDim, fieldVtsOffsets[ifld]);
    ifld++;
  }

  message += indent(3) + "</PointData>\n";
  message += indent(2) + "</Piece>\n";
  message += indent(1) + "</StructuredGrid>\n";
  message += indent(1) + "<AppendedData encoding=\"raw\">\n";
  message += "_";

  const MPI_Offset appendedStart = static_cast<MPI_Offset>(message.size());

  if (rank == 0) {
    MPI_File_write_at(file_out, 0, message.c_str(), message.size(), MPI_CHAR, MPI_STATUS_IGNORE);
  }

  MPI_Barrier(mpi_comm);

  // ---------------------------------------------------------------------------
  // Helper: write one appended VTK block.
  //
  // vtkOffset is relative to appendedStart, i.e. immediately after "_".
  // ---------------------------------------------------------------------------

  auto writeField = [&](const MPI_Offset vtkOffset, const int nComponents, const std::vector<float> &field) {
    const unsigned long long int nbytes = static_cast<unsigned long long int>(nComponents) *
                                          static_cast<unsigned long long int>(v.nGlobal) * sizeof(float);

    const MPI_Offset blockStart = appendedStart + vtkOffset;
    const MPI_Offset dataStart = blockStart + static_cast<MPI_Offset>(sizeof(unsigned long long int));

    const MPI_Offset localDataStart = dataStart + static_cast<MPI_Offset>(sizeof(float)) *
                                                      static_cast<MPI_Offset>(nComponents) *
                                                      static_cast<MPI_Offset>(pOffset);

    if (rank == 0) {
      MPI_File_write_at(file_out, blockStart, &nbytes, 1, MPI_UNSIGNED_LONG_LONG, MPI_STATUS_IGNORE);
    }

    MPI_Barrier(mpi_comm);

    MPI_File_write_at_all(file_out, localDataStart, field.data(), field.size(), MPI_FLOAT, MPI_STATUS_IGNORE);
  };

  // Write coordinates.
  {
    std::vector<float> positions(dim * v.nLocal, 0.0f);

    for (dlong i = 0; i < v.nLocal; i++) {
      positions[dim * i + 0] = static_cast<float>(v.x[i]);
      positions[dim * i + 1] = static_cast<float>(v.y[i]);
      positions[dim * i + 2] = static_cast<float>(v.z[i]);
    }

    writeField(posVtsOffset, dim, positions);
  }

  // Write one appended block per VTK DataArray.
  ifld = 0;

  for (const auto &entry : userFieldList) {
    const auto &fieldVector = std::get<1>(entry);
    const int fieldDim = fieldVector.size();

    std::vector<float> oneField(fieldDim * v.nLocal);

    for (int d = 0; d < fieldDim; d++) {
      for (dlong i = 0; i < v.nLocal; i++) {
        oneField[d * v.nLocal + i] = v.data[v.offsetScan[ifld] + d * v.nLocal + i];
      }
    }

    writeField(fieldVtsOffsets[ifld], fieldDim, oneField);

    ifld++;
  }

  // Write XML footer at the exact end of appended data.
  const std::string footer = "\n  </AppendedData>\n"
                             "</VTKFile>\n";

  const MPI_Offset footerOffset = appendedStart + static_cast<MPI_Offset>(totalAppendedBytes);

  if (rank == 0) {
    MPI_File_write_at(file_out, footerOffset, footer.c_str(), footer.size(), MPI_CHAR, MPI_STATUS_IGNORE);
  }

  MPI_Barrier(mpi_comm);

  MPI_File_close(&file_out);
}

inline void fieldExtract::process(const double time, const int tstep)
{
  nekrsCheck(!setupCalled, MPI_COMM_SELF, EXIT_FAILURE, "%s\n", "process called prior to setup!");

  interpolate();

  stepCounter++;
  const int out_step = stepCounter;

  const int nx = boxDims[0];
  const int ny = (boxDim > 1) ? boxDims[1] : 1;
  const int nz = (boxDim > 2) ? boxDims[2] : 1;
  const std::string tag = (boxDim == 1) ? "Line" : (boxDim == 2) ? "Plane" : "Box";

  GridView v{std::array<int, 3>{nx, ny, nz},
             numPoints,
             numPointsLocal,
             numPointsScan,
             xCoord,
             yCoord,
             zCoord,
             fldData,
             fldDataOffsetScan,
             point0,
             point1};

  std::ostringstream output;
  output << fileName << tag << std::setw(5) << std::setfill('0') << out_step << ".vts";
  writeVts(v, output.str(), time, tstep);
}

inline void
fieldExtract::doAvg(const double time, const int tstep, const std::string &avgDir, const std::string &mode)
{
  nekrsCheck(!setupCalled, MPI_COMM_SELF, EXIT_FAILURE, "%s\n", "doAvg called prior to setup!");
  nekrsCheck(!boxMode, MPI_COMM_SELF, EXIT_FAILURE, "%s\n", "doAvg requires box mode (no custom points)!");
  nekrsCheck(mode != "gather" && mode != "gather-scatter",
             MPI_COMM_SELF,
             EXIT_FAILURE,
             "invalid doAvg mode '%s' (use gather|gather-scatter)\n",
             mode.c_str());

  // parse avgDir -> which axes to average
  std::array<bool, 3> avg{false, false, false};
  {
    int bad = 0;
    for (char c : avgDir) {
      const int a = (c == 'x' || c == 'X') ? 0 : (c == 'y' || c == 'Y') ? 1 : (c == 'z' || c == 'Z') ? 2 : -1;
      if (a < 0 || avg[a]) {
        bad = 1;
        break;
      }
      avg[a] = true;
    }
    const bool any = avg[0] || avg[1] || avg[2];
    nekrsCheck(bad || !any,
               MPI_COMM_SELF,
               EXIT_FAILURE,
               "invalid avgDir '%s' (non-repeating subset of x,y,z)\n",
               avgDir.c_str());
  }

  const int n[3] = {boxDims[0], (boxDim > 1) ? boxDims[1] : 1, (boxDim > 2) ? boxDims[2] : 1};
  for (int a = 0; a < 3; a++) {
    nekrsCheck(avg[a] && n[a] < 2,
               MPI_COMM_SELF,
               EXIT_FAILURE,
               "cannot average singleton axis '%c'\n",
               "xyz"[a]);
  }

  const int nx = n[0], ny = n[1], nz = n[2];
  const int nxy = nx * ny;

  // canonical direction tag (x,y,z order)
  std::string dirTag;
  for (int a = 0; a < 3; a++) {
    if (avg[a]) {
      dirTag += "xyz"[a];
    }
  }

  // reduced grid: averaged axes collapse to 1
  const int rn[3] = {avg[0] ? 1 : nx, avg[1] ? 1 : ny, avg[2] ? 1 : nz};
  const dlong rN = static_cast<dlong>(rn[0]) * rn[1] * rn[2];
  const int rnxy = rn[0] * rn[1];

  auto &ags = setupAvgGs(avg, dirTag); // cached after the first call per direction
  const dlong nTotal = ags.nTotal;
  const dlong rLocal = ags.rLocal;
  const dlong rScan = ags.rScan;

  // fresh field values on the box
  interpolate();

  int nComp = 0; // total number of scalar components over all fields
  for (auto &entry : userFieldList) {
    nComp += std::get<1>(entry).size();
  }

  auto trap = [](int i, int nn) { return (nn > 1 && (i == 0 || i == nn - 1)) ? 0.5 : 1.0; };

  // local pre-reduction: accumulate trapezoid-weighted partial sums into the
  // per-cell slots (one gs vector of length nTotal per component); the
  // reduced-grid tail entries contribute zero
  std::vector<dfloat> wrk(static_cast<size_t>(nComp) * nTotal, 0.0);
  {
    int comp = 0, ifld = 0;
    for (const auto &entry : userFieldList) {
      const int fieldDim = std::get<1>(entry).size();
      for (int d = 0; d < fieldDim; d++) {
        for (dlong i = 0; i < numPointsLocal; i++) {
          const dlong ig = i + numPointsScan;
          int ix, iy, iz;
          getIxyz(ig, nx, nxy, ix, iy, iz);
          dfloat w = 1.0;
          if (avg[0]) {
            w *= trap(ix, nx);
          }
          if (avg[1]) {
            w *= trap(iy, ny);
          }
          if (avg[2]) {
            w *= trap(iz, nz);
          }
          wrk[static_cast<size_t>(comp) * nTotal + ags.pointSlot[i]] +=
              w * fldData[fldDataOffsetScan[ifld] + d * numPointsLocal + i];
        }
        comp++;
      }
      ifld++;
    }
  }

  // sum the partial sums of all ranks sharing a reduced cell; the total is
  // scattered back to every member (slots AND reduced-grid tail)
  oogs::startFinish(wrk.data(), nComp, nTotal, ogsDfloat, ogsAdd, ags.gsh);

  // trapezoidal normalization: product of (n_a - 1) over averaged axes
  dfloat denom = 1.0;
  for (int a = 0; a < 3; a++) {
    if (avg[a]) {
      denom *= (n[a] - 1);
    }
  }
  for (auto &s : wrk) {
    s /= denom;
  }

  // output stream name + per-stream counter
  const std::string stem = fileName + "_avg" + dirTag;
  const int out_step = ++avgCounter[stem + "_" + mode];

  if (mode == "gather-scatter") {
    // scatter each point's cell average back through pointSlot
    // (reuse box coords/dims/split)
    std::vector<float> sdata(fldDataOffsetScan[Nfields]);
    int comp = 0, ifld = 0;
    for (const auto &entry : userFieldList) {
      const int fieldDim = std::get<1>(entry).size();
      for (int d = 0; d < fieldDim; d++) {
        for (dlong i = 0; i < numPointsLocal; i++) {
          sdata[fldDataOffsetScan[ifld] + d * numPointsLocal + i] =
              static_cast<float>(wrk[static_cast<size_t>(comp) * nTotal + ags.pointSlot[i]]);
        }
        comp++;
      }
      ifld++;
    }

    const std::string tag = (boxDim == 1) ? "Line" : (boxDim == 2) ? "Plane" : "Box";
    GridView v{std::array<int, 3>{nx, ny, nz},
               numPoints,
               numPointsLocal,
               numPointsScan,
               xCoord,
               yCoord,
               zCoord,
               sdata,
               fldDataOffsetScan,
               point0,
               point1};
    std::ostringstream output;
    output << stem << tag << std::setw(5) << std::setfill('0') << out_step << ".vts";
    writeVts(v, output.str(), time, tstep);
    return;
  }

  // mode == "gather": emit the collapsed grid; the tail entries of wrk are the
  // reduced grid, already in the pointDistribution(rN) split across ranks
  const dfloat mid[3] = {static_cast<dfloat>(0.5) * (point0[0] + point1[0]),
                         static_cast<dfloat>(0.5) * (point0[1] + point1[1]),
                         static_cast<dfloat>(0.5) * (point0[2] + point1[2])};
  const dfloat dcoord[3] = {(nx > 1) ? (point1[0] - point0[0]) / (nx - 1) : 0.0,
                            (ny > 1) ? (point1[1] - point0[1]) / (ny - 1) : 0.0,
                            (nz > 1) ? (point1[2] - point0[2]) / (nz - 1) : 0.0};

  std::vector<dfloat> rx(rLocal), ry(rLocal), rz(rLocal);
  for (dlong i = 0; i < rLocal; i++) {
    const dlong ig = i + rScan;
    int jx, jy, jz;
    getIxyz(ig, rn[0], rnxy, jx, jy, jz);
    rx[i] = avg[0] ? mid[0] : point0[0] + jx * dcoord[0];
    ry[i] = avg[1] ? mid[1] : point0[1] + jy * dcoord[1];
    rz[i] = avg[2] ? mid[2] : point0[2] + jz * dcoord[2];
  }

  // local reduced data slice + its offset scan
  std::vector<size_t> rLocalOffset(Nfields + 1);
  rLocalOffset[0] = 0;
  {
    int ifld = 1;
    for (auto &entry : userFieldList) {
      const auto dimf = std::get<1>(entry).size();
      rLocalOffset[ifld] = rLocalOffset[ifld - 1] + static_cast<size_t>(rLocal) * dimf;
      ifld++;
    }
  }
  std::vector<float> rdata(rLocalOffset[Nfields]);
  {
    int comp = 0, ifld = 0;
    for (const auto &entry : userFieldList) {
      const int fieldDim = std::get<1>(entry).size();
      for (int d = 0; d < fieldDim; d++) {
        for (dlong i = 0; i < rLocal; i++) {
          rdata[rLocalOffset[ifld] + d * rLocal + i] =
              static_cast<float>(wrk[static_cast<size_t>(comp) * nTotal + ags.nSlots + i]);
        }
        comp++;
      }
      ifld++;
    }
  }

  int eff = 0;
  for (int a = 0; a < 3; a++) {
    if (rn[a] > 1) {
      eff++;
    }
  }
  const std::string tag = (eff == 0) ? "Point" : (eff == 1) ? "Line" : (eff == 2) ? "Plane" : "Box";

  GridView v{std::array<int, 3>{rn[0], rn[1], rn[2]},
             rN,
             rLocal,
             rScan,
             rx,
             ry,
             rz,
             rdata,
             rLocalOffset,
             point0,
             point1};
  std::ostringstream output;
  output << stem << tag << std::setw(5) << std::setfill('0') << out_step << ".vts";
  writeVts(v, output.str(), time, tstep);
}

#endif // nekrs_field_extract_hpp_
