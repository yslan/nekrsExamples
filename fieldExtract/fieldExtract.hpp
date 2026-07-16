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
// Grid points follow lexicographic (nz, ny, nx) ordering and are distributed
// evenly across MPI ranks. The .vts layout is consumed by test.py.
//
// Roadmap (see plans/): cylinder / Gauss-Radau layouts, derived fields (dT/dz),
// redistribution of mismatched explicit points, in-memory data() accessor.

#if !defined(nekrs_field_extract_hpp_)
#define nekrs_field_extract_hpp_

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

  ~fieldExtract() = default;

  void process(const double time, const int tstep); // interpolate + write .vts

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

  bool setupCalled = false;
  bool setPointsCalled = false;

  void getIxyz(const dlong i, const int nx, const int nxy, int &ix, int &iy, int &iz);
  void setBoxDims(const std::vector<int> &boxDims_in);

  void setupCommon(const std::vector<int> &boxDims_in,
                   const std::vector<field> &fldList,
                   const std::string &fileName_in);
  void finishSetup(mesh_t *mesh);

  void setPoints(const std::array<dfloat, 3> &x0, const std::array<dfloat, 3> &x1);
  void setPoints(const std::array<std::vector<dfloat>, 3> &XYZ);

  void interpolate();
  void writeFld(const std::string fname, const double time, const int tstep);
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
  setupCommon(boxDims_in, fldList, fileName_in);
  // point0/point1 stay {0,0,0}; they only feed the .vts x0/x1 tags.
  setPoints(XYZ);
  finishSetup(mesh);
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

inline void fieldExtract::writeFld(const std::string fname, const double time, const int tstep)
{
  MPI_Comm mpi_comm = platform->comm.mpiComm();
  const int rank = platform->comm.mpiRank();

  // Local point offset in the global point ordering.
  const long long int pOffset = numPointsScan;

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

  const int nx = (boxDim >= 1) ? boxDims[0] : 1;
  const int ny = (boxDim >= 2) ? boxDims[1] : 1;
  const int nz = (boxDim >= 3) ? boxDims[2] : 1;
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
  vtsOffset += headerBytes + dim * numPoints * sizeof(float);

  std::vector<long long int> fieldVtsOffsets;
  fieldVtsOffsets.reserve(Nfields);

  for (const auto &entry : userFieldList) {
    const auto &fieldVector = std::get<1>(entry);
    const int fieldDim = fieldVector.size();

    fieldVtsOffsets.push_back(vtsOffset);
    vtsOffset += headerBytes + fieldDim * numPoints * sizeof(float);
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
  message += indent(3) + fieldDataArrayInt64("numPoints", static_cast<long long>(numPoints));
  message += indent(3) + fieldDataArrayFloat3("x0", point0[0], point0[1], point0[2]);
  message += indent(3) + fieldDataArrayFloat3("x1", point1[0], point1[1], point1[2]);
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
                                          static_cast<unsigned long long int>(numPoints) * sizeof(float);

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
    std::vector<float> positions(dim * numPointsLocal, 0.0f);

    for (dlong i = 0; i < numPointsLocal; i++) {
      positions[dim * i + 0] = static_cast<float>(xCoord[i]);
      positions[dim * i + 1] = static_cast<float>(yCoord[i]);
      positions[dim * i + 2] = static_cast<float>(zCoord[i]);
    }

    writeField(posVtsOffset, dim, positions);
  }

  // Write one appended block per VTK DataArray.
  ifld = 0;

  for (const auto &entry : userFieldList) {
    const auto &fieldVector = std::get<1>(entry);
    const int fieldDim = fieldVector.size();

    std::vector<float> oneField(fieldDim * numPointsLocal);

    for (int d = 0; d < fieldDim; d++) {
      for (dlong i = 0; i < numPointsLocal; i++) {
        oneField[d * numPointsLocal + i] = fldData[fldDataOffsetScan[ifld] + d * numPointsLocal + i];
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

  { // write vtk
    std::string tag;
    if (boxDim == 1) {
      tag = "Line";
    } else if (boxDim == 2) {
      tag = "Plane";
    } else {
      tag = "Box";
    }

    std::ostringstream output;
    output << fileName << tag << std::setw(5) << std::setfill('0') << out_step << ".vts";
    std::string fname = output.str();

    writeFld(fname, time, tstep);
  }
}

#endif // nekrs_field_extract_hpp_
