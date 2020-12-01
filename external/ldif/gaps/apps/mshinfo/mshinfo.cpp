// Source file for the mesh info program



// Include files 

namespace gaps {}
using namespace gaps;
#include "R3Shapes/R3Shapes.h"



int NumBoundaryEdges(R3Mesh *mesh)
{
  // Count number of boundary edges
  int count = 0;
  for (int i = 0; i < mesh->NEdges(); i++) {
    R3MeshEdge *edge = mesh->Edge(i);
    R3MeshFace *face0 = mesh->FaceOnEdge(edge, 0);
    R3MeshFace *face1 = mesh->FaceOnEdge(edge, 1);
    if (!face0 || !face1) count++;
  }
  return count;
}



int main(int argc, char **argv)
{
  // Check number of arguments
  if (argc != 2) {
    printf("Usage: meshinfo inputfile\n");
    exit(1);
  }

  // Create mesh
  R3Mesh *mesh = new R3Mesh();
  if (!mesh) {
    RNFail("Unable to allocate mesh data structure.\n");
    exit(1);
  }

  // Read mesh
  if (!mesh->ReadFile(argv[1])) {
    RNFail("Unable to read file: %s\n", argv[1]);
    exit(1);
  }

  // Compute stats
  int num_boundaries = NumBoundaryEdges(mesh);
  int num_components = mesh->ConnectedComponents();
  RNScalar avg_vertex_valence = mesh->AverageVertexValence();
  RNScalar avg_edge_length = mesh->AverageEdgeLength();
  RNArea total_face_area = mesh->Area();
  RNAngle avg_edge_interior_angle = mesh->AverageEdgeInteriorAngle();
  const R3Point& centroid = mesh->Centroid();
  const R3Box& bbox = mesh->BBox();

  // Print stats
  printf("Total faces = %d\n", mesh->NFaces());
  printf("Total edges = %d\n", mesh->NEdges());
  printf("Total vertices = %d\n", mesh->NVertices());
  printf("Boundary edges = %d ( %g%% )\n", num_boundaries, 100.0 * num_boundaries/ mesh->NEdges());
  printf("Connected components = %d\n", num_components);
  printf("Surface area = %g\n", total_face_area);
  printf("Average vertex valence = %g\n", avg_vertex_valence);
  printf("Average edge length = %g\n", avg_edge_length);
  printf("Average edge interior angle = %g\n", avg_edge_interior_angle);
  printf("Centroid = ( %g %g %g )\n", centroid[0], centroid[1], centroid[2]);
  printf("Bounding box = ( %g %g %g ) ( %g %g %g )\n", bbox[0][0], bbox[0][1], bbox[0][2], bbox[1][0], bbox[1][1], bbox[1][2]);
  printf("Axial lengths = ( %g %g %g )\n", bbox.XLength(), bbox.YLength(), bbox.ZLength());

  // Return success 
  return 0;
}

















