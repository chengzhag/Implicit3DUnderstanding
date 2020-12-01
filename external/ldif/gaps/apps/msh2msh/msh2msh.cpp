// Source file for the mesh converter program



// Include files 

namespace gaps {}
using namespace gaps;
#include "R3Shapes/R3Shapes.h"



// Program arguments

const char *input_mesh_name = NULL;
const char *output_mesh_name = NULL;
const char *output_xform_name = NULL;
const char *color_name = NULL;
const char *merge_list_name = NULL;
int flip_faces = 0;
int clean = 0;
int swap_edges = 0;
int fill_holes = 0;
int delete_interior_faces = 0;
RNScalar smooth_factor = 0;
R3Affine xform(R4Matrix(1,0,0,0, 0,1,0,0, 0,0,1,0, 0,0,0,1));
RNLength min_edge_length = 0;
RNLength max_edge_length = 0;
RNLength min_component_area = 0;
char *xform_name = NULL;
int scale_by_area = 0;
int scale_by_pca = 0;
int rotate_by_pca = 0;
int align_by_pca = 0;
int translate_by_centroid = 0;
int print_verbose = 0;



////////////////////////////////////////////////////////////////////////
// I/O STUFF
////////////////////////////////////////////////////////////////////////

static R3Mesh *
ReadMesh(const char *mesh_name)
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

  // Check if mesh is valid
  assert(mesh->IsValid());

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



static int
WriteMesh(R3Mesh *mesh, const char *filename)
{
  // Start statistics
  RNTime start_time;
  start_time.Read();

  // Write mesh to file
  if (!mesh->WriteFile(filename)) return 0;

  // Print statistics
  if (print_verbose) {
    printf("Wrote mesh to %s ...\n", filename);
    printf("  Time = %.2f seconds\n", start_time.Elapsed());
    printf("  # Faces = %d\n", mesh->NFaces());
    printf("  # Edges = %d\n", mesh->NEdges());
    printf("  # Vertices = %d\n", mesh->NVertices());
    fflush(stdout);
  }

  // Return success
  return 1;
}



static int
ReadMatrix(R4Matrix& m, const char *filename)
{
  // Open file
  FILE *fp = fopen(filename, "r");
  if (!fp) {
    RNFail("Unable to open matrix file: %s\n", filename);
    return 0;
  }

  // Read matrix from file
  for (int i = 0; i < 4; i++) {
    for (int j = 0; j < 4; j++) {
      double value;
      fscanf(fp, "%lf", &value);
      m[i][j] = value;
    }
  }

  // Close file
  fclose(fp);

  // Return success
  return 1;
}



static int
WriteMatrix(const R4Matrix& m, const char *filename)
{
  // Open file
  FILE *fp = fopen(filename, "w");
  if (!fp) {
    RNFail("Unable to open matrix file: %s\n", filename);
    return 0;
  }

  // Write matrix to file
  for (int i = 0; i < 4; i++) {
    for (int j = 0; j < 4; j++) {
      fprintf(fp, "%g ", m[i][j]);
    }
  }

  // Close file
  fclose(fp);

  // Return success
  return 1;
}



////////////////////////////////////////////////////////////////////////
// PROCESSING STUFF
////////////////////////////////////////////////////////////////////////

static int
MergeList(R3Mesh *mesh, const char *filename)
{
  // Open list file
  FILE *fp = fopen(filename, "r");
  if (!fp) {
    RNFail("Unable to open list %s\n", filename);
    return 0;
  }

  // Merge meshes from list
  char src_filename[4096];
  while (fscanf(fp, "%s", src_filename) == (unsigned int) 1) {
    // Read mesh
    R3Mesh src_mesh;
    if (!src_mesh.ReadFile(src_filename)) return 0;
    mesh->CreateCopy(src_mesh);
  }

  // Close file
  fclose(fp);

  // Return success
  return 1;
}



static int
CopyColors(R3Mesh *mesh, const char *source_mesh_name)
{
  // Read source mesh
  R3Mesh source_mesh;
  if (!source_mesh.ReadFile(source_mesh_name)) return 0;

  // Create kdtree
  R3MeshSearchTree kdtree(&source_mesh);
  
  // Copy colors
  for (int i = 0; i < mesh->NVertices(); i++) {
    R3MeshVertex *vertex = mesh->Vertex(i);
    mesh->SetVertexColor(vertex, RNblack_rgb);
    const R3Point& position = mesh->VertexPosition(vertex);

    // Search kdtree
    R3MeshIntersection closest;
    kdtree.FindClosest(position, closest);
    if (closest.type == R3_MESH_VERTEX_TYPE) {
      mesh->SetVertexColor(vertex, source_mesh.VertexColor(closest.vertex));
    }
    else if (closest.type == R3_MESH_EDGE_TYPE) {
      R3Span span = mesh->EdgeSpan(closest.edge);
      RNScalar t = span.T(position);
      int k = (t < 0.5 * span.Length()) ?  0 : 1;
      R3MeshVertex *source_vertex = source_mesh.VertexOnEdge(closest.edge, k);
      mesh->SetVertexColor(vertex, source_mesh.VertexColor(source_vertex));
    }
    else if (closest.type == R3_MESH_FACE_TYPE) {
      R3Point b = source_mesh.FaceBarycentric(closest.face, position);
      int k = b.Vector().MaxDimension();
      R3MeshVertex *source_vertex = source_mesh.VertexOnFace(closest.face, k);
      mesh->SetVertexColor(vertex, source_mesh.VertexColor(source_vertex));
    }
  }

  // Return success
  return 1;
}



static int
DeleteSmallConnectedComponents(R3Mesh *mesh, RNArea min_area)
{
  // Start statistics
  RNTime start_time;
  start_time.Read();
  int total_ncomponents = 0;
  int small_ncomponents = 0;

  // Allocate seed marks
  if (mesh->NFaces() == 0) return 0;
  int *components = new int [ mesh->NFaces() ];
  for (int i = 0; i < mesh->NFaces(); i++) components[i]  = -1;

  // Iterate finding connected components
  int seed_index = 0;
  RNArray<R3MeshFace *> small_component_faces;
  for (;;) {
    // Find an unmarked face
    R3MeshFace *seed = NULL;
    for (; seed_index < mesh->NFaces(); seed_index++) {
      R3MeshFace *face = mesh->Face(seed_index);
      if (components[seed_index] < 0) {
        seed = face;
        break;
      }
    }

    // Check if found a new component
    if (!seed) break;

    // Mark connected component
    R3mesh_mark++;
    RNArea area = 0;
    RNArray<R3MeshFace *> stack;
    RNArray<R3MeshFace *> component;
    stack.InsertTail(seed);
    mesh->SetFaceMark(seed, R3mesh_mark);
    while (!stack.IsEmpty()) {
      R3MeshFace *face = stack.Tail();
      stack.RemoveTail();

      // Add face to component
      component.Insert(face);
      area += mesh->FaceArea(face);
      components[mesh->FaceID(face)] = total_ncomponents;

      // Check neighbors
      for (int i = 0; i < 3; i++) {
        R3MeshFace *neighbor = mesh->FaceOnFace(face, i);
        if (!neighbor) continue;
        if (mesh->FaceMark(neighbor) == R3mesh_mark) continue;
        mesh->SetFaceMark(neighbor, R3mesh_mark);
        stack.InsertTail(neighbor);
      }
    }

    // Check component
    if (component.NEntries() == 0) continue;

    // Mark component for deletion if too small
    if (area < min_area) {
      small_component_faces.Append(component);
      small_ncomponents++;
    }

    // Increment number of components
    total_ncomponents++;
  }

  // Delete deletable faces
  for (int i = 0; i < small_component_faces.NEntries(); i++) {
    R3MeshFace *face = small_component_faces.Kth(i);
    mesh->DeleteFace(face);
  }

  // Delete components
  delete [] components;

  // Print debug statistics
  if (print_verbose) {
    printf("  Deleted small connected components ...\n");
    printf("    Time = %.2f seconds\n", start_time.Elapsed());
    printf("    # Deleted Components = %d\n", small_ncomponents);
    printf("    # Deleted Faces = %d\n", small_component_faces.NEntries());
    printf("    # Remaining Components = %d\n", total_ncomponents - small_ncomponents);
    printf("    # Remaining Faces = %d\n", mesh->NFaces());
    fflush(stdout);
  }

  // Return success
  return 1;
}



static int
DeleteInteriorFaces(R3Mesh *mesh)
{
  // Start statistics
  RNTime start_time;
  start_time.Read();

  // Get convenient variables
  int resolution = 1024;
  R3Point c = mesh->BBox().Centroid();
  RNScalar r = mesh->BBox().DiagonalRadius();
  R2Box bbox2d(-r, -r, r, r);
  R2Grid face_image(bbox2d, r / resolution);
  R2Grid depth_image(face_image);
  RNScalar d = face_image.GridToWorldScaleFactor();
  RNArea pixel_area = d*d;
  const int nviews = 100;

  // Initialize face marks
  int *marks = new int [ mesh->NFaces() ];
  for (int i = 0; i < mesh->NFaces(); i++) marks[i] = 0;

  // Mark faces visible from any exterior view
  for (int k = 0; k < nviews; k++) {
    // Create transformation
    R3Affine transformation = R3identity_affine;
    transformation.Rotate(R3RandomDirection(), RN_TWO_PI * RNRandomScalar());
    transformation.Rotate(R3RandomDirection(), RN_TWO_PI * RNRandomScalar());
    transformation.Rotate(R3RandomDirection(), RN_TWO_PI * RNRandomScalar());
    transformation.Translate(-c.Vector());

    // Render image
    face_image.Clear(-1);
    depth_image.Clear(RN_INFINITY);
    for (int i = 0; i < mesh->NFaces(); i++) {
      R3MeshFace *face = mesh->Face(i);
      R3MeshVertex *v0 = mesh->VertexOnFace(face, 0);
      R3MeshVertex *v1 = mesh->VertexOnFace(face, 1);
      R3MeshVertex *v2 = mesh->VertexOnFace(face, 2);
      R3Point p0 = mesh->VertexPosition(v0);
      R3Point p1 = mesh->VertexPosition(v1);
      R3Point p2 = mesh->VertexPosition(v2);
      p0.Transform(transformation);
      p1.Transform(transformation);
      p2.Transform(transformation);
      face_image.RenderWorldTriangle(
        R2Point(p0.X(), p0.Y()),
        R2Point(p1.X(), p1.Y()),
        R2Point(p2.X(), p2.Y()),
        i, i, i,
        p0.Z(), p1.Z(), p2.Z(),
        depth_image);
    }

    // Mark visible faces
    for (int i = 0; i < face_image.NEntries(); i++) {
      RNScalar face_value = face_image.GridValue(i);
      if (face_value < 0) continue;
      int face_index = (int) (face_value + 0.5);
      if (face_index < 0) continue;
      if (face_index >= mesh->NFaces()) continue;
      marks[face_index]++;
    }
  }

  // Find unmarked faces
  RNArray<R3MeshFace *> unmarked_faces;
  for (int i = 0; i < mesh->NFaces(); i++) {
    R3MeshFace *face = mesh->Face(i);
    int possible_marks = nviews * mesh->FaceArea(face) / pixel_area;
    if (marks[i] > 1E-3*possible_marks) continue;
    unmarked_faces.Insert(face);
  }
  
  // Delete unmarked faces
  for (int i = 0; i < unmarked_faces.NEntries(); i++) {
    R3MeshFace *face = unmarked_faces.Kth(i);
    mesh->DeleteFace(face);
  }

  // Print debug statistics
  if (print_verbose) {
    printf("  Deleted interior faces ...\n");
    printf("    Time = %.2f seconds\n", start_time.Elapsed());
    printf("    # Deleted Faces = %d\n", unmarked_faces.NEntries());
    printf("    # Remaining Faces = %d\n", mesh->NFaces());
    fflush(stdout);
  }

  // Return success
  return 1;
}



////////////////////////////////////////////////////////////////////////
// PROGRAM ARGUMENT PARSING
////////////////////////////////////////////////////////////////////////

int ParseArgs(int argc, char **argv)
{
  // Check number of arguments
  if (argc == 1) {
    printf("Usage: mesh2mesh inputname outputname [options]\n");
    exit(0);
  }

  // Parse arguments
  argc--; argv++;
  while (argc > 0) {
    if ((*argv)[0] == '-') {
      R3Affine prev_xform = xform;
      if (!strcmp(*argv, "-v")) print_verbose = 1;
      else if (!strcmp(*argv, "-flip")) flip_faces = 1;
      else if (!strcmp(*argv, "-clean")) clean = 1;
      else if (!strcmp(*argv, "-fill_holes")) fill_holes = 1;
      else if (!strcmp(*argv, "-delete_interior_faces")) delete_interior_faces = 1;
      else if (!strcmp(*argv, "-swap_edges")) swap_edges = 1;
      else if (!strcmp(*argv, "-scale_by_area")) scale_by_area = 1;
      else if (!strcmp(*argv, "-scale_by_pca")) scale_by_pca = 1;
      else if (!strcmp(*argv, "-translate_by_centroid")) translate_by_centroid = 1;
      else if (!strcmp(*argv, "-center_at_origin")) translate_by_centroid = 1;
      else if (!strcmp(*argv, "-rotate_by_pca")) rotate_by_pca = 1;
      else if (!strcmp(*argv, "-align_by_pca")) align_by_pca = 1;
      else if (!strcmp(*argv, "-smooth"))  { argv++; argc--; smooth_factor = atof(*argv); }
      else if (!strcmp(*argv, "-scale")) { argv++; argc--; xform = R3identity_affine; xform.Scale(atof(*argv)); xform.Transform(prev_xform); }
      else if (!strcmp(*argv, "-tx")) { argv++; argc--; xform = R3identity_affine; xform.XTranslate(atof(*argv)); xform.Transform(prev_xform); }
      else if (!strcmp(*argv, "-ty")) { argv++; argc--; xform = R3identity_affine; xform.YTranslate(atof(*argv)); xform.Transform(prev_xform);}
      else if (!strcmp(*argv, "-tz")) { argv++; argc--; xform = R3identity_affine; xform.ZTranslate(atof(*argv)); xform.Transform(prev_xform);}
      else if (!strcmp(*argv, "-sx")) { argv++; argc--; xform = R3identity_affine; xform.XScale(atof(*argv)); xform.Transform(prev_xform);}
      else if (!strcmp(*argv, "-sy")) { argv++; argc--; xform = R3identity_affine; xform.YScale(atof(*argv)); xform.Transform(prev_xform);}
      else if (!strcmp(*argv, "-sz")) { argv++; argc--; xform = R3identity_affine; xform.ZScale(atof(*argv)); xform.Transform(prev_xform);}
      else if (!strcmp(*argv, "-rx")) { argv++; argc--; xform = R3identity_affine; xform.XRotate(RN_PI*atof(*argv)/180.0); xform.Transform(prev_xform);}
      else if (!strcmp(*argv, "-ry")) { argv++; argc--; xform = R3identity_affine; xform.YRotate(RN_PI*atof(*argv)/180.0); xform.Transform(prev_xform);}
      else if (!strcmp(*argv, "-rz")) { argv++; argc--; xform = R3identity_affine; xform.ZRotate(RN_PI*atof(*argv)/180.0); xform.Transform(prev_xform);}
      else if (!strcmp(*argv, "-xform")) { argv++; argc--; R4Matrix m;  if (ReadMatrix(m, *argv)) { xform = R3identity_affine; xform.Transform(R3Affine(m)); xform.Transform(prev_xform);} } 
      else if (!strcmp(*argv, "-min_edge_length")) { argv++; argc--; min_edge_length = atof(*argv); }
      else if (!strcmp(*argv, "-max_edge_length")) { argv++; argc--; max_edge_length = atof(*argv); }
      else if (!strcmp(*argv, "-remove_small_components")) { argv++; argc--; min_component_area = atof(*argv); }
      else if (!strcmp(*argv, "-color")) { argv++; argc--; color_name = *argv; }
      else if (!strcmp(*argv, "-merge_list")) { argv++; argc--; merge_list_name = *argv; }
      else if (!strcmp(*argv, "-debug_matrix")) { argv++; argc--; output_xform_name = *argv; }
      else if (!strcmp(*argv, "-transform")) {
        R4Matrix m;
        argv++; argc--; m[0][0] = atof(*argv); argv++; argc--; m[0][1] = atof(*argv);
        argv++; argc--; m[0][2] = atof(*argv); argv++; argc--; m[0][3] = atof(*argv);
        argv++; argc--; m[1][0] = atof(*argv); argv++; argc--; m[1][1] = atof(*argv);
        argv++; argc--; m[1][2] = atof(*argv); argv++; argc--; m[1][3] = atof(*argv);
        argv++; argc--; m[2][0] = atof(*argv); argv++; argc--; m[2][1] = atof(*argv);
        argv++; argc--; m[2][2] = atof(*argv); argv++; argc--; m[2][3] = atof(*argv);
        argv++; argc--; m[3][0] = atof(*argv); argv++; argc--; m[3][1] = atof(*argv);
        argv++; argc--; m[3][2] = atof(*argv); argv++; argc--; m[3][3] = atof(*argv);
        xform = R3identity_affine;
        xform.Transform(R3Affine(m));
        xform.Transform(prev_xform);
      }  
      else { RNFail("Invalid program argument: %s", *argv); exit(1); }
      argv++; argc--;
    }
    else {
      if (!input_mesh_name) input_mesh_name = *argv;
      else if (!output_mesh_name) output_mesh_name = *argv;
      else { RNFail("Invalid program argument: %s", *argv); exit(1); }
      argv++; argc--;
    }
  }

  // Check input filename
  if (!input_mesh_name) {
    RNFail("You did not specify an input file name.\n");
    return 0;
  }

  // Check output filename
  if (!output_mesh_name) {
    RNFail("You did not specify an output file name.\n");
    return 0;
  }

  // Return OK status 
  return 1;
}



////////////////////////////////////////////////////////////////////////
// MAIN
////////////////////////////////////////////////////////////////////////

int main(int argc, char **argv)
{
  // Check number of arguments
  if (!ParseArgs(argc, argv)) exit(1);

  // Read mesh
  R3Mesh *mesh = ReadMesh(input_mesh_name);
  if (!mesh) exit(-1);

  // Merge meshes from list
  if (merge_list_name) {
    if (!MergeList(mesh, merge_list_name)) exit(-1);
  }

  // Rotate so that principal axes are aligned with canonical coordinate frame
  if (rotate_by_pca) {
    R3Point centroid = mesh->Centroid();
    R3Triad triad = mesh->PrincipleAxes(&centroid);
    R3Affine rotation = R3identity_affine;
    rotation.Translate(centroid.Vector());
    rotation.Transform(R3Affine(triad.InverseMatrix()));
    rotation.Translate(-centroid.Vector());
    xform.Transform(rotation);
  }

  // Scale based on mesh area
  if (scale_by_area) {
    RNArea area = mesh->Area();
    if (area > 0) xform.Scale(1 / sqrt(area));
  }

  // Scale based on principal component analysis
  if (scale_by_pca) {
    RNScalar variances[3];
    R3Point centroid = mesh->Centroid();
    mesh->PrincipleAxes(&centroid, variances);
    if (variances[0] > 0) xform.Scale(1.0 / sqrt(variances[0]));
  }

  // Center mesh at origin
  if (translate_by_centroid) {
    // R3Point centroid = R3zero_point;
    // for (int i = 0; i < mesh->NVertices(); i++) {
    //   R3MeshVertex *vertex = mesh->Vertex(i);
    //   centroid += mesh->VertexPosition(vertex);
    // }
    // if (mesh->NVertices() > 0) centroid /= mesh->NVertices();
    R3Point centroid = mesh->Centroid();
    xform.Translate(-centroid.Vector());
  }

  // Normalize translation, rotation, and scale
  if (align_by_pca) {
    xform = mesh->PCANormalizationTransformation();
  }

  // Transform 
  if (!xform.IsIdentity()) {
    mesh->Transform(xform);
  }

  // Delete interior faces
  if (delete_interior_faces) {
    DeleteInteriorFaces(mesh);
  }

  // Remove small components
  if (min_component_area > 0) {
    DeleteSmallConnectedComponents(mesh, min_component_area);
  }

  // Clean 
  if (clean) {
    mesh->DeleteUnusedEdges();
    mesh->DeleteUnusedVertices();
  }

  // Flip every face
  if (flip_faces) {
    for (int i = 0; i < mesh->NFaces(); i++) {
      mesh->FlipFace(mesh->Face(i));
    }
  }

  // Fill holes
  if (fill_holes) {
    mesh->FillHoles();
  }

  // Smooth
  if (smooth_factor > 0) {
    mesh->Smooth(smooth_factor);
  }

  // Subdivide edges that are too long
  if (max_edge_length > 0) {
    mesh->SubdivideLongEdges(max_edge_length);
  }

  // Split edges that are too long
  if (min_edge_length > 0) {
    mesh->CollapseShortEdges(min_edge_length);
  }

  // Swap edges
  if (swap_edges) {
    mesh->SwapEdges();
  }

  // Transfer colors
  if (color_name) {
    CopyColors(mesh, color_name);
  }
  
  // Write mesh
  if (!WriteMesh(mesh, output_mesh_name)) exit(-1);

  // Write transformation
  if (!xform.IsIdentity() && output_xform_name) {
    if (!WriteMatrix(xform.Matrix(), output_xform_name)) exit(-1);
  }

  // Return success 
  return 0;
}




