/*
 * Copyright (C) 2015   Malte Brunn
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

//------------------------------------------------------------------------------
#include "typedef.hpp"
#include "communicator.hpp"
#include "compute.hpp"
#include "geometry.hpp"
#include "parameter.hpp"
#include "grid.hpp"

#ifdef USE_DEBUG_VISU
#include "visu.hpp"
#endif // USE_DEBUG_VISU
#ifdef USE_VTK
#include "vtk.hpp"
#include <sys/stat.h>
#endif // USE_VTK

#include "argvparser.hpp"
#include "precice/SolverInterface.hpp"

using namespace precice;
using namespace precice::constants;

int main(int argc, char **argv) {
  Communicator comm(&argc, &argv);
  Parameter param;
  Geometry geom(&comm);
  ARGVParser parser;
  parser.bind("-geom", [&geom](int ac, char **av) -> int {
    if (ac != 1) return 0;
    geom.Load(av[0]);
    return 1;
  });
  parser.bind("-param", [&param](int ac, char **av) -> int {
    if (ac != 1) return 0;
    param.Load(av[0]);
    return 1;
  });
  parser.exec(argc, argv);
  // Create the fluid solver
  Compute comp(&geom, &param);
  // Create parameter and geometry instances with default values
  int N = geom.Coup(); // Number of mesh elements
  if(!param.Expl()){
      N *= 4;
  }
  std::string config("./precice-config.xml");
  std::string solverName("Fluid");
  SolverInterface interface(solverName,0,1);

  interface.configure(config);
  // get domain dimension from precice
  int dim = interface.getDimensions();
  std::string mesh_name("Fluid-Mesh");
  int meshID = interface.getMeshID(mesh_name);
  int temperatureID = interface.getDataID("Temperature", meshID);
  int heatfluxID = interface.getDataID("Heat-Flux", meshID);
  int *vertexIDs = new int[N];
  double *vertices = new double[dim * N];
  for(int i=0;i<N;i++){
      vertices[i]=-1;
  }
  double *temperature, *heatflux;
  temperature = new double[N];
  heatflux = new double[N];
  comp.Vertices(vertices,temperature,heatflux,N,dim);

  interface.setMeshVertices(meshID, N, vertices, vertexIDs);

  std::cout << "Initialize preCICE..." << std::endl;
  double precice_dt = interface.initialize();
  const std::string& coric = actionReadIterationCheckpoint();
  const std::string& cowic = actionWriteIterationCheckpoint();
  if (interface.isActionRequired(actionWriteInitialData())) {
      interface.writeBlockScalarData(temperatureID,N,vertexIDs,temperature); // write initial Temperature
      interface.fulfilledAction(actionWriteInitialData());
  }
  interface.initializeData(); // synchronize with OpenFOAM
  if (interface.isReadDataAvailable()) {
      interface.readBlockScalarData(heatfluxID,N,vertexIDs,heatflux); // read heatfluxCoupled
  }

#ifdef USE_VTK
  if (comm.getRank() == 0) {
    // check if folder "VTK" exists
    struct stat info;

    if (stat("VTK", &info) != 0) {
      system("mkdir VTK");
    }
  }
#endif

// Create and initialize the visualization
#ifdef USE_DEBUG_VISU
  Renderer visu(geom.Length(), geom.Mesh());
  double ratio = geom.Length()[1]/geom.Length()[0];
  visu.Init(800/ comm.ThreadDim()[0], 800*ratio/ comm.ThreadDim()[1], comm.getRank() + 1);
#endif // USE_DEBUG_VISU

multi_real_t offset;
offset[0] = comm.ThreadIdx()[0] * (geom.Mesh()[0] * (double)(geom.Size()[0] - 2));
offset[1] = comm.ThreadIdx()[1] * (geom.Mesh()[1] * (double)(geom.Size()[1] - 2));
#ifdef USE_VTK
VTK vtk(geom.Mesh(), geom.Length(), geom.TotalLength(), offset, comm.getRank(),
      comm.getSize(), comm.ThreadDim());
#endif

#ifdef USE_DEBUG_VISU
  const Grid *visugrid;

  visugrid = comp.GetVelocity();
#endif // USE_DEBUG_VISU

while (interface.isCouplingOngoing()) {
    #ifdef USE_DEBUG_VISU
    // Render and check if window is closed
    switch (visu.Render(visugrid)) {
    case -1:
      return -1;
    case 0:
      visugrid = comp.GetVelocity();
      break;
    case 1:
      visugrid = comp.GetU();
      break;
    case 2:
      visugrid = comp.GetV();
      break;
    case 3:
      visugrid = comp.GetP();
      break;
    case 4:
      visugrid = comp.GetT();
    default:
      break;
    };
    #endif // USE_DEBUG_VISU

    #ifdef USE_VTK
    // Create VTK Files in the folder VTK
    vtk.Init("VTK/field");
    vtk.AddRank();
    vtk.AddCellField("Cell Velocity", comp.GetU(), comp.GetV());
    vtk.SwitchToPointData();
    vtk.AddPointField("Velocity", comp.GetU(), comp.GetV());
    vtk.AddPointScalar("Pressure", comp.GetP());
    vtk.AddPointScalar("Temperature", comp.GetT());
    vtk.Finish();
    #endif
    for(int iCount=0; iCount <1; iCount++){
        comp.TimeStep(true,&interface, temperatureID, heatfluxID,N,vertexIDs,
    vertices,temperature,heatflux,precice_dt,coric,cowic);
    }

  }
  interface.finalize();
  return 0;
}
