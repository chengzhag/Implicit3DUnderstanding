// Source file for the mesh to grid conversion program



// Include files 

namespace gaps {}
using namespace gaps;
#include "R3Shapes/R3Shapes.h"



// Program variables

static char *mesh_name = NULL;
static char *grid_name = NULL;
static int min_resolution = 16;
static int max_resolution = 512;
static double grid_spacing = 0.1;
static int sdf = 0;
static int print_verbose = 0;



static R3Mesh *
ReadMesh(char *mesh_name)
{
  // Start statistics
  RNTime start_time;
  start_time.Read();

  // Allocate mesh
  R3Mesh *mesh = new R3Mesh();
  assert(mesh);

  // Read mesh from file
  if (!mesh->ReadFile(mesh_name)) {
    delete mesh;
    return NULL;
  }

  // Print statistics
  if (print_verbose) {
    printf("Read mesh ...\n");
    printf("  Time = %.2f seconds\n", start_time.Elapsed());
    printf("  # Faces = %d\n", mesh->NFaces());
    printf("  # Edges = %d\n", mesh->NEdges());
    printf("  # Vertices = %d\n", mesh->NVertices());
    fflush(stdout);
  }

  // Return success
  return mesh;
}



static RNLength
ShortestEdgeLength(R3Mesh *mesh, R3MeshFace *face)
{
  // Return shortest edge on face
  RNLength d0 = mesh->EdgeLength(mesh->EdgeOnFace(face, 0));
  RNLength d1 = mesh->EdgeLength(mesh->EdgeOnFace(face, 1));
  RNLength d2 = mesh->EdgeLength(mesh->EdgeOnFace(face, 2));
  RNLength d = FLT_MAX;
  if (d0 < d) d = d0;
  if (d1 < d) d = d1;
  if (d2 < d) d = d2;
  return d;
}



static R3Grid *
CreateGrid(R3Mesh *mesh)
{
  // Start statistics
  RNTime start_time;
  start_time.Read();

  // Allocate grid
  R3Grid *grid = new R3Grid(mesh->BBox(), grid_spacing, min_resolution, max_resolution);
  if (!grid) {
    RNFail("Unable to allocate grid\n");
    exit(-1);
  }

  // Rasterize each triangle into grid
  for (int i = 0; i < mesh->NFaces(); i++) {
    R3MeshFace *face = mesh->Face(i);
    const R3Point& p0 = mesh->VertexPosition(mesh->VertexOnFace(face, 0));
    const R3Point& p1 = mesh->VertexPosition(mesh->VertexOnFace(face, 1));
    const R3Point& p2 = mesh->VertexPosition(mesh->VertexOnFace(face, 2));
    grid->RasterizeWorldTriangle(p0, p1, p2, 1.0);
  }

  // Make binary (in case triangles overlap)
  grid->Threshold(0.5, 0, 1);

  // Check if should compute sdf
  if (sdf) {
    // Compute distance grid
    grid->SquaredDistanceTransform();
    grid->Sqrt();
    grid->Multiply(grid->GridToWorldScaleFactor());

    // Allocate temmporary memory
    int *votes = new int [ grid->NEntries() ];
    int *components = new int [ grid->NEntries() ];
    for (int i = 0; i < grid->NEntries(); i++) votes[i] = 0;
    for (int i = 0; i < grid->NEntries(); i++) components[i] = 0;

    // Vote for signs by connected component
    R3Grid signs(*grid);
    signs.Clear(0);
    for (int i = 0; i < mesh->NFaces(); i++) {
      R3MeshFace *face = mesh->Face(i);
      RNLength e = ShortestEdgeLength(mesh, face);
      if (e > grid_spacing) e = grid_spacing;
      R3Point p = mesh->FaceCentroid(face);
      R3Point p1 = p + e * mesh->FaceNormal(face);
      R3Point p2 = p - e * mesh->FaceNormal(face);
      signs.RasterizeWorldPoint(p1, 1.0);
      signs.RasterizeWorldPoint(p2, -1.0);
    }

    // Vote for connected components
    int ncomponents = grid->ConnectedComponents(0, grid->NEntries(), NULL, NULL, components);
    if (ncomponents > 0) {
      for (int i = 0; i < grid->NEntries(); i++) {
        int component = components[i];
        if (component < 0) continue;
        votes[component] += signs.GridValue(i);
      }

      // Set signs by votes on connected components
      for (int i = 0; i < grid->NEntries(); i++) {
        int component = components[i];
        if (component < 0) continue;
        if (votes[component] < 0) {
          RNScalar d = grid->GridValue(i);
          grid->SetGridValue(i, -d);
        }
      }
    }
    
    // Delete temporary memory
    delete [] components;
    delete [] votes;
  }

  // Print statistics
  if (print_verbose) {
    printf("Created grid ...\n");
    printf("  Time = %.2f seconds\n", start_time.Elapsed());
    printf("  Resolution = %d %d %d\n", grid->XResolution(), grid->YResolution(), grid->ZResolution());
    printf("  Spacing = %g\n", grid->GridToWorldScaleFactor());
    printf("  Cardinality = %d\n", grid->Cardinality());
    printf("  Volume = %g\n", grid->Volume());
    RNInterval grid_range = grid->Range();
    printf("  Minimum = %g\n", grid_range.Min());
    printf("  Maximum = %g\n", grid_range.Max());
    printf("  L1Norm = %g\n", grid->L1Norm());
    printf("  L2Norm = %g\n", grid->L2Norm());
    fflush(stdout);
  }

  // Return grid
  return grid;
}



static int 
WriteGrid(R3Grid *grid, const char *grid_name)
{
  // Start statistics
  RNTime start_time;
  start_time.Read();

  // Write grid
  int status = grid->WriteFile(grid_name);

  // Print statistics
  if (print_verbose) {
    printf("Wrote grid ...\n");
    printf("  Time = %.2f seconds\n", start_time.Elapsed());
    printf("  # Bytes = %d\n", status * (int) sizeof(RNScalar));
    fflush(stdout);
  }

  // Return status
  return status;
}



static int 
ParseArgs(int argc, char **argv)
{
  // Parse arguments
  argc--; argv++;
  while (argc > 0) {
    if ((*argv)[0] == '-') {
      if (!strcmp(*argv, "-v")) print_verbose = 1; 
      else if (!strcmp(*argv, "-sdf")) sdf = 1; 
      else if (!strcmp(*argv, "-spacing")) { argc--; argv++; grid_spacing = atof(*argv); }
      else if (!strcmp(*argv, "-min_resolution")) { argc--; argv++; min_resolution = atoi(*argv); }
      else if (!strcmp(*argv, "-max_resolution")) { argc--; argv++; max_resolution = atoi(*argv); }
      else { RNFail("Invalid program argument: %s", *argv); exit(1); }
    }
    else {
      if (!mesh_name) mesh_name = *argv;
      else if (!grid_name) grid_name = *argv;
      else { RNFail("Invalid program argument: %s", *argv); exit(1); }
    }
    argv++; argc--;
  }

  // Check mesh filename
  if (!mesh_name || !grid_name) {
    RNFail("Usage: msh2grd mesh grid [options]\n");
    return 0;
  }

  // Return OK status 
  return 1;
}



int 
main(int argc, char **argv)
{
  // Parse program arguments
  if (!ParseArgs(argc, argv)) exit(-1);

  // Read mesh file
  R3Mesh *mesh = ReadMesh(mesh_name);
  if (!mesh) exit(-1);

  // Create grid from mesh
  R3Grid *grid = CreateGrid(mesh);
  if (!grid) exit(-1);

  // Write grid
  int status = WriteGrid(grid, grid_name);
  if (!status) exit(-1);

  // Return success
  return 0;
}
