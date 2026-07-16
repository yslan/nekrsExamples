#if !defined(nekrs_field_extract_hpp_)
#define nekrs_field_extract_hpp_

#include <variant>
#include <functional>
#include <optional>

#include "platform.hpp"
#include "mesh3D.h"
#include "pointInterpolation.hpp"

class fieldExtract
{
public:

  using field = std::tuple<std::string, std::vector<occa::memory>>;

  fieldExtract(mesh_t *mesh_in,
               const std::vector<int> &boxDims_in,
               const std::vector<field> &_in,
               const std::string fileName_in, // TODO: maybe put into ->process?, tavg->writeToFile for IO
               const std::array<dfloat, 3>& x0 = std::array<dfloat, 3>{0,0,0}, // TODO: std::vector< tuple<min,max> >, check size
               const std::array<dfloat, 3>& x1 = std::array<dfloat, 3>{0,0,0}, 
               const std::array<std::vector<dfloat>, 3>& XYZ = std::array<std::vector<dfloat>, 3>{});

  // cylinder
  //Gauss-Radau in r, uniform coordinates

  ~fieldExtract(); // TODO:

//  std::vector<dfloat> data(); // TODO

  void process(const double time, const int tstep); // interpolation, seperate writeToFile

private:
  pointInterpolation_t *interpolator = nullptr;
//  comm_t comm; // TODO neknek?

  std::vector<field> userFieldList;
  std::vector<size_t> fieldOffsetScan;
  int Nfields;

  std::string fileName;
  int stepCounter;

  int boxDim = 0;
  std::vector<int> boxDims;

  std::array<dfloat, 3> point0; // start pt
  std::array<dfloat, 3> point1; // end pt

  dlong numPoints;
  dlong numPointsLocal;
  dlong numPointsScan;

  std::vector<dfloat> xCoord, yCoord, zCoord;
  dlong pointsOffset;
  std::vector<size_t> pointsOffsetScan;

  std::vector<float> fldData;


//TODO save coordinates? see vtk writer
  bool initialized = false;
  bool setupCalled = false;
  bool setPointsCalled = false;

  void initialize();
  void get_ixyz(const int i, const int nx, const int nxy, int &ix, int &iy, int &iz);
  void set_boxDims(const std::vector<int> &boxDims_in);
  void set_connectivity();

  void setPoints(const std::array<dfloat, 3>& x0, const std::array<dfloat, 3>& x1);
  void setPoints(const std::array<std::vector<dfloat>, 3>& XYZ);

  void interpolate();
  void writeFld(const std::string fname, const double time, const int tstep);
};

#endif // hpp

void fieldExtract::initialize()
{
  platform->timer.addUserStat("fieldExtract::"); 
  initialized = true;// TODO
}

void fieldExtract::get_ixyz(const int i, const int nx, const int nxy, int &ix, int &iy, int &iz)
{
  iz = i / nxy;
  int rem = i % nxy;
  iy = rem / nx;
  ix = rem % nx;
}

void fieldExtract::set_boxDims(const std::vector<int> &boxDims_in)
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
  nekrsCheck(numPoints < 1,
             MPI_COMM_SELF,
             EXIT_FAILURE,
             "invalid number of points: %d!\n",
             numPoints);
  
  const int rank = platform->comm.mpiRank();
  const int Nrank = platform->comm.mpiCommSize();
  const int npts_base = numPoints / Nrank;   // first evenly distribute points
  const int npts_left = numPoints % Nrank;   // assign remaining points
  numPointsLocal = (rank < npts_left) ? npts_base + 1 : npts_base;
  
  numPointsScan = 0;
  if (rank > 0) { 
    numPointsScan = rank * npts_base + std::min(npts_left, rank);
  }
}

void fieldExtract::setPoints(const std::array<dfloat, 3>& point0, const std::array<dfloat, 3>& point1)
{
  // setup points from x0 to x1
  if (numPointsLocal) {
    const int nx = boxDims[0];
    const int ny = (boxDim > 1) ? boxDims[1] : 1;
    const int nz = (boxDim > 2) ? boxDims[2] : 1;
    const int nxy = nx * ny;

    const dfloat dx = (nx > 1) ? (point1[0] - point0[0]) / (nx - 1) : 0.0;
    const dfloat dy = (ny > 1) ? (point1[1] - point0[1]) / (ny - 1) : 0.0;
    const dfloat dz = (nz > 1) ? (point1[2] - point0[2]) / (nz - 1) : 0.0;

    xCoord.resize(numPointsLocal);
    yCoord.resize(numPointsLocal);
    zCoord.resize(numPointsLocal);

    for (int i=0; i<numPointsLocal; i++) {
      const int ig = i + numPointsScan;
      int ix, iy, iz;
      get_ixyz(ig, nx, nxy, ix, iy, iz);
      xCoord[i] = point0[0] + ix * dx;
      yCoord[i] = point0[1] + iy * dy;
      zCoord[i] = point0[2] + iz * dz;
    }
  }

  setPointsCalled = true;
}

void fieldExtract::setPoints(const std::array<std::vector<dfloat>, 3>& XYZ)
{

  dlong Nlocal = XYZ[0].size();

  { // chk XYZ
    int err = 0;
    for (int d = 1; d < 3; d++) {
      if (XYZ[d].size() != Nlocal) {
        err = d;
        break;
      }
    }

    MPI_Allreduce(MPI_IN_PLACE, &err, 1, MPI_INT, MPI_MAX, platform->comm.mpiComm()); // for err-print

    nekrsCheck(err != 0,
               MPI_COMM_SELF,
               EXIT_FAILURE,
               "XYZ[%d] has different length as XYZ[0]!\n",
               err);
  }

  { // chk XYZ vs boxDims
    nekrsCheck(Nlocal != numPointsLocal, // TODO: redistribute points
               MPI_COMM_SELF,
               EXIT_FAILURE,
               "input points don't match the boxDims! Ntotal=%d  numPoints=%d\n",
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

  for (int i = 0; i < numPointsLocal; ++i) {
    xCoord[i] = XYZ[0][i];
    yCoord[i] = XYZ[1][i];
    zCoord[i] = XYZ[2][i];
  }

  setPointsCalled = true;
}


fieldExtract::fieldExtract(mesh_t *mesh_in,
                           const std::vector<int> &boxDims_in,
                           const std::vector<field> &userFieldList_in,
                           const std::string fileName_in,
                           const std::array<dfloat, 3>& point0_in,
                           const std::array<dfloat, 3>& point1_in,
                           const std::array<std::vector<dfloat>, 3> &XYZ)
{
  if (!initialized) {
    initialize();
  }

  nekrsCheck(setupCalled,
             MPI_COMM_SELF,
             EXIT_FAILURE,
             "%s\n",
             "Attempt to call setup twice!\n");

  auto mesh = mesh_in;

  set_boxDims(boxDims_in);

  userFieldList = userFieldList_in;

  {
    Nfields = userFieldList.size();
    fieldOffsetScan.resize(Nfields + 1);

    int ifld = 1; 
    fieldOffsetScan[0] = 0;
    for (auto &entry : userFieldList) {
      auto o_fld = std::get<1>(entry);
      const auto dim_fld = o_fld.size();

      auto len_fld = o_fld[0].size();

      // len_fld must be the same for all dim
      int err = 0;
      for (int idim = 1; idim < dim_fld; idim++) {
        if (len_fld != o_fld[0].size()) {
          err = idim;
          break;
        }
      }

      MPI_Allreduce(MPI_IN_PLACE, &err, 1, MPI_INT, MPI_MAX, platform->comm.mpiComm()); // for err-print
      nekrsCheck(err != 0,
                 MPI_COMM_SELF,
                 EXIT_FAILURE,
                 "userFieldList[%d][%d] has different length than its 0th array!\n",
                 ifld, err);

      // TODO: should we read fieldOffset?
      fieldOffsetScan[ifld] =
          fieldOffsetScan[ifld - 1] + alignStride<dfloat>(len_fld) * dim_fld;

      ifld++;  
    }
  }

  fileName = fileName_in;
  stepCounter = 0;

  // setup points
  point0 = point0_in;
  point1 = point1_in;
  int mode1 = 0;
  {
    dfloat dist = 0;
    for (int d = 0; d < mesh->dim; d++) {
      dist += pow(point0[d] - point1[d], 2);
    }
    if (dist > 1e-6) mode1 = 1;
    MPI_Allreduce(MPI_IN_PLACE, &mode1, 1, MPI_INT, MPI_MAX, platform->comm.mpiComm());
  }

  int mode2 = 1;
  {
    for (int d = 0; d < mesh->dim; d++) {
      if (XYZ[d].empty()) {
        mode2 = 0;  
        break;
      }
    }
    MPI_Allreduce(MPI_IN_PLACE, &mode2, 1, MPI_INT, MPI_MAX, platform->comm.mpiComm());    
  }

  nekrsCheck(mode1 > 0 && mode2 > 0,
             MPI_COMM_SELF,
             EXIT_FAILURE,
             "%s\n",
             "Cannot specify both endpoints and points!");

  if (mode1 > 0) {
    setPoints(point0, point1);
  }
  if (mode2 > 0) {
    setPoints(XYZ);
  }

  nekrsCheck(!setPointsCalled,
             MPI_COMM_SELF,
             EXIT_FAILURE,
             "%s\n", 
             "points are not setup properly!");

  pointsOffset = numPointsLocal;
  pointsOffsetScan.resize(Nfields + 1);

  int ifld = 1; 
  pointsOffsetScan[0] = 0;
  for (auto &entry : userFieldList) {
    auto o_fld = std::get<1>(entry);
    const auto dim_fld = o_fld.size();
    
    // TODO: should we read fieldOffset?
    pointsOffsetScan[ifld] =
        pointsOffsetScan[ifld - 1] + pointsOffset * dim_fld;
    
    ifld++;
  }

  fldData.resize(pointsOffsetScan[Nfields]);

  interpolator = new pointInterpolation_t(mesh, platform->comm.mpiComm());
  interpolator->setPoints(xCoord, yCoord, zCoord);
  interpolator->find();

  setupCalled = true;
}

void fieldExtract::interpolate()
{
  platform->timer.tic("fieldExtractnek::interpolate");
  auto o_wrk = platform->deviceMemoryPool.reserve<dfloat>(pointsOffset);
  std::vector<dfloat> wrk(pointsOffset);

  int ifld = 0;
  for (const auto& entry : userFieldList) {
    auto fieldVector = std::get<1>(entry);
    auto fieldDim = fieldVector.size();

    for (int d = 0; d < fieldDim; d++) {
      interpolator->eval(1, fieldOffsetScan[ifld], fieldVector[d], pointsOffset, o_wrk); // todo vec?
      o_wrk.copyTo(wrk.data(), pointsOffset);

      for (dlong i = 0; i < numPointsLocal; i++) {
        fldData[pointsOffsetScan[ifld] + d * pointsOffset + i] = static_cast<float>(wrk[i]);
      }
    }
    ifld++;
  }
  platform->timer.toc("fieldExtractnek::interpolate");
}

void fieldExtract::writeFld(const std::string fname,
                            const double time,
                            const int tstep)
{
  MPI_Comm mpi_comm = platform->comm.mpiComm();
  const int rank = platform->comm.mpiRank();

  // Local point offset in the global point ordering.
  long long int pOffset = numPointsScan;
/*  MPI_Exscan(&numPointsLocal,
             &pOffset,
             1,
             MPI_LONG_LONG_INT,
             MPI_SUM,
             mpi_comm);
  if (rank == 0) {
    pOffset = 0;
  }*/

  // Remove old file / truncate cleanly before MPI-IO opens it.
  if (rank == 0) {
    std::ofstream file(fname, std::ios::binary | std::ios::trunc);
    file.close();
  }
  MPI_Barrier(mpi_comm);

  MPI_File file_out;
  MPI_File_open(mpi_comm,
                fname.c_str(),
                MPI_MODE_CREATE | MPI_MODE_WRONLY,
                MPI_INFO_NULL,
                &file_out);

  auto indent = [](const int n)
  {
    return std::string(2 * n, ' ');
  };

  auto fieldDataArrayInt3 = [](const std::string& name,
                               const int a,
                               const int b,
                               const int c)
  {
    return "<DataArray type=\"Int32\" Name=\"" + name +
           "\" NumberOfComponents=\"3\" NumberOfTuples=\"1\" format=\"ascii\"> " +
           std::to_string(a) + " " +
           std::to_string(b) + " " +
           std::to_string(c) +
           " </DataArray>\n";
  };
  
  auto fieldDataArrayFloat3 = [](const std::string& name,
                                 const double a,
                                 const double b,
                                 const double c)
  {
    return "<DataArray type=\"Float64\" Name=\"" + name +
           "\" NumberOfComponents=\"3\" NumberOfTuples=\"1\" format=\"ascii\"> " +
           std::to_string(a) + " " +
           std::to_string(b) + " " +
           std::to_string(c) +
           " </DataArray>\n";
  };
  
  auto fieldDataArrayInt64 = [](const std::string& name,
                                const long long v)
  {
    return "<DataArray type=\"Int64\" Name=\"" + name +
           "\" NumberOfTuples=\"1\" format=\"ascii\"> " +
           std::to_string(v) +
           " </DataArray>\n";
  };

  auto vtsDataArray = [](const std::string &fieldName,
                         long long int nComponent,
                         long long int offset)
  {
    return "<DataArray type=\"Float32\" Name=\"" + fieldName +
           "\" NumberOfComponents=\"" + std::to_string(nComponent) +
           "\" format=\"appended\" offset=\"" +
           std::to_string(offset) + "\"/>\n";
  };

  const int nx = (boxDim >= 1) ? boxDims[0] : 1;
  const int ny = (boxDim >= 2) ? boxDims[1] : 1;
  const int nz = (boxDim >= 3) ? boxDims[2] : 1;
  const std::string vtsExtent =
      "0 " + std::to_string(nx - 1) + " " +
      "0 " + std::to_string(ny - 1) + " " +
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

  for (const auto& entry : userFieldList) {
    const auto& fieldVector = std::get<1>(entry);
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
  message += indent(3) + "<DataArray type=\"Float32\" Name=\"TimeValue\" "
                         "NumberOfTuples=\"1\" format=\"ascii\"> " +
                         std::to_string(time) + " </DataArray>\n";
  message += indent(3) + "<DataArray type=\"Int32\" Name=\"CYCLE\" "
                         "NumberOfTuples=\"1\" format=\"ascii\"> " +
                         std::to_string(tstep) + " </DataArray>\n";
  message += indent(3) + fieldDataArrayInt3("boxDims", nx, ny, nz);
  message += indent(3) +
             fieldDataArrayInt64("numPoints",
                                 static_cast<long long>(numPoints));
  message += indent(3) +
             fieldDataArrayFloat3("x0", point0[0], point0[1], point0[2]);
  message += indent(3) +
             fieldDataArrayFloat3("x1", point1[0], point1[1], point1[2]);
  message += indent(2) + "</FieldData>\n";

  message += indent(2) + "<Piece Extent=\"" + vtsExtent + "\">\n";
  message += indent(3) + "<Points>\n";
  message += indent(4) + vtsDataArray("Position", dim, posVtsOffset);
  message += indent(3) + "</Points>\n";
  message += indent(3) + "<PointData>\n";

  int ifld = 0;
  for (const auto& entry : userFieldList) {
    const std::string fieldName = std::get<0>(entry);
    const auto& fieldVector = std::get<1>(entry);
    const int fieldDim = fieldVector.size();
    message += indent(4) +
               vtsDataArray(fieldName, fieldDim, fieldVtsOffsets[ifld]);
    ifld++;
  }

  message += indent(3) + "</PointData>\n";
  message += indent(2) + "</Piece>\n";
  message += indent(1) + "</StructuredGrid>\n";
  message += indent(1) + "<AppendedData encoding=\"raw\">\n";
  message += "_";

  const MPI_Offset appendedStart =
      static_cast<MPI_Offset>(message.size());

  if (rank == 0) {
    MPI_File_write_at(file_out,
                      0,
                      message.c_str(),
                      message.size(),
                      MPI_CHAR,
                      MPI_STATUS_IGNORE);
  }

  MPI_Barrier(mpi_comm);

  // ---------------------------------------------------------------------------
  // Helper: write one appended VTK block.
  //
  // vtkOffset is relative to appendedStart, i.e. immediately after "_".
  // ---------------------------------------------------------------------------

  auto writeField = [&](const MPI_Offset vtkOffset,
                        const int nComponents,
                        const std::vector<float>& field)
  {
    const unsigned long long int nbytes =
        static_cast<unsigned long long int>(nComponents) *
        static_cast<unsigned long long int>(numPoints) *
        sizeof(float);

    const MPI_Offset blockStart = appendedStart + vtkOffset;
    const MPI_Offset dataStart =
        blockStart + static_cast<MPI_Offset>(sizeof(unsigned long long int));

    const MPI_Offset localDataStart =
        dataStart
        + static_cast<MPI_Offset>(sizeof(float))
        * static_cast<MPI_Offset>(nComponents)
        * static_cast<MPI_Offset>(pOffset);

    if (rank == 0) {
      MPI_File_write_at(file_out,
                        blockStart,
                        &nbytes,
                        1,
                        MPI_UNSIGNED_LONG_LONG,
                        MPI_STATUS_IGNORE);
    }

    MPI_Barrier(mpi_comm);

    MPI_File_write_at_all(file_out,
                          localDataStart,
                          field.data(),
                          field.size(),
                          MPI_FLOAT,
                          MPI_STATUS_IGNORE);
  };

  // Write coordinates.
  {
    std::vector<float> positions(dim * pointsOffset, 0.0f);

    for (dlong i = 0; i < numPointsLocal; i++) {
      positions[dim * i + 0] = static_cast<float>(xCoord[i]);
      positions[dim * i + 1] = static_cast<float>(yCoord[i]);
      positions[dim * i + 2] = static_cast<float>(zCoord[i]);
    }

    writeField(posVtsOffset, dim, positions);
  }

  // Write one appended block per VTK DataArray.
  ifld = 0;

  for (const auto& entry : userFieldList) {
    const auto& fieldVector = std::get<1>(entry);
    const int fieldDim = fieldVector.size();

    std::vector<float> oneField(fieldDim * pointsOffset);

    for (int d = 0; d < fieldDim; d++) {
      for (dlong i = 0; i < pointsOffset; i++) {
        oneField[d * pointsOffset + i] =
            fldData[pointsOffsetScan[ifld] + d * pointsOffset + i];
      }
    }

    writeField(fieldVtsOffsets[ifld], fieldDim, oneField);

    ifld++;
  }

  // Write XML footer at the exact end of appended data.
  const std::string footer =
      "\n  </AppendedData>\n"
      "</VTKFile>\n";

  const MPI_Offset footerOffset =
      appendedStart + static_cast<MPI_Offset>(totalAppendedBytes);

  if (rank == 0) {
    MPI_File_write_at(file_out,
                      footerOffset,
                      footer.c_str(),
                      footer.size(),
                      MPI_CHAR,
                      MPI_STATUS_IGNORE);
  }

  MPI_Barrier(mpi_comm);

  MPI_File_close(&file_out);
}

void fieldExtract::process(const double time, const int tstep)
{
  nekrsCheck(!setupCalled,
             MPI_COMM_SELF,
             EXIT_FAILURE,
             "%s\n",
             "process called prior to setup!");

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

/*
std::vector<float> fieldExtract::data()
{
  return ;
}
*/

fieldExtract::~fieldExtract()
{
/*
  userFieldList.clear();//TODO
  o_AVG.free();
  fldWriter.reset();*/
}
