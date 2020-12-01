// Source file for the mesh segmentation program

////////////////////////////////////////////////////////////////////////
// Include files
////////////////////////////////////////////////////////////////////////

namespace gaps {}
using namespace gaps;
#include "R3Shapes/R3Shapes.h"



////////////////////////////////////////////////////////////////////////
// Program arguments
////////////////////////////////////////////////////////////////////////

static char *input_mesh_filename = NULL;
static char *output_mesh_filename = NULL;
static char *output_texture_directory = NULL;
static int separate_face_segments = 0;
static int separate_face_materials = 0;
static double texel_spacing = 1;
static int flatten = 0;
static int print_verbose = 0;



////////////////////////////////////////////////////////////////////////
// Type definitions
////////////////////////////////////////////////////////////////////////

struct Segment {
  R3Mesh *mesh;
  RNArray<R3MeshFace *> faces;
  RNArray<R3MeshVertex *> vertices;
  RNMap<R3MeshVertex *, R3MeshVertex *> vertex_map;
  R4Matrix matrix;
  R2Box bbox;
  int id;
  int index;
};



////////////////////////////////////////////////////////////////////////
// Input/output stuff
////////////////////////////////////////////////////////////////////////

static R3Mesh *
ReadMesh(const char *filename)
{
  // Start statistics
  RNTime start_time;
  start_time.Read();

  // Allocate mesh
  R3Mesh *mesh = new R3Mesh();
  assert(mesh);

  // Read mesh from file
  if (!mesh->ReadFile(filename)) {
    delete mesh;
    return NULL;
  }

  // Print statistics
  if (print_verbose) {
    printf("Read mesh from %s ...\n", filename);
    printf("  Time = %.2f seconds\n", start_time.Elapsed());
    printf("  # Faces = %d\n", mesh->NFaces());
    printf("  # Edges = %d\n", mesh->NEdges());
    printf("  # Vertices = %d\n", mesh->NVertices());
    fflush(stdout);
  }

  // Return mesh
  return mesh;
}



static int
WriteMesh(R3Mesh *mesh, char *filename)
{
  // Start statistics
  RNTime start_time;
  start_time.Read();

  // Write mesh
  if (!mesh->WriteFile(filename)) {
    fprintf(stderr, "Unable to write mesh to %s\n", filename);
    return 0;
  }

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



////////////////////////////////////////////////////////////////////////
// Segmentation stuff
////////////////////////////////////////////////////////////////////////

static int 
CompareSegments(const void *data1, const void *data2)
{
  // Sort by bbox.xlength
  Segment *segment1 = *((Segment **) data1);
  Segment *segment2 = *((Segment **) data2);
  if (segment1->bbox.XLength() > segment2->bbox.XLength()) return -1;
  else if (segment1->bbox.XLength() < segment2->bbox.XLength()) return 1;
  else return 0;
}
  


static int 
CreateSegments(R3Mesh *mesh, RNArray<Segment *>& segments)
{
  // Remember original vertices
  RNArray<R3MeshVertex *> original_vertices;
  for (int i = 0; i < mesh->NVertices(); i++) {
    R3MeshVertex *original_vertex = mesh->Vertex(i);
    original_vertices.Insert(original_vertex);
  }

  // Remember original faces
  RNArray<R3MeshFace *> original_faces;
  for (int i = 0; i < mesh->NFaces(); i++) {
    R3MeshFace *original_face = mesh->Face(i);
    original_faces.Insert(original_face);
  }
  
  // Create segments
  RNMap<int, Segment *> segment_map;
  R3MeshVertex *segment_face_vertices[3];
  for (int i = 0; i < original_faces.NEntries(); i++) {
    R3MeshFace *original_face = original_faces.Kth(i);
    int id = 0;
    if (separate_face_segments) id = mesh->FaceSegment(original_face);
    else if (separate_face_materials) id = mesh->FaceMaterial(original_face);

    // Find segment
    Segment *segment = NULL;
    if (!segment_map.Find(id, &segment)) {
      segment = new Segment();
      segment->mesh = mesh;
      segment->matrix = R4identity_matrix;
      segment->bbox = R2null_box;
      segment->id = id;
      segment->index = segments.NEntries();
      segment_map.Insert(id, segment);
      segments.Insert(segment);
    }

    // Find/insert copy of vertices into segment
    for (int j = 0; j < 3; j++) {
      R3MeshVertex *original_vertex = mesh->VertexOnFace(original_face, j);

      // Find/insert vertex in segment
      R3MeshVertex *segment_vertex = NULL;
      if (!segment->vertex_map.Find(original_vertex, &segment_vertex)) {
        segment_vertex = mesh->CreateVertex(*original_vertex);
        assert(segment_vertex);
        segment->vertex_map.Insert(original_vertex, segment_vertex);
        segment->vertices.Insert(segment_vertex);
      }

      // Remember vertex for face
      segment_face_vertices[j] = segment_vertex;
    }
      
    // Insert copy of face into segment
    R3MeshFace *segment_face = mesh->CreateFace(segment_face_vertices[0], segment_face_vertices[1], segment_face_vertices[2]);
    assert(segment_face);
    mesh->SetFaceCategory(segment_face, mesh->FaceCategory(original_face));
    mesh->SetFaceSegment(segment_face, mesh->FaceSegment(original_face));
    mesh->SetFaceMaterial(segment_face, mesh->FaceMaterial(original_face));
    segment->faces.Insert(segment_face);
  }

  // Empty segment vertex maps (because will delete original vertices next)
  for (int i = 0; i < segments.NEntries(); i++) {
    Segment *segment = segments.Kth(i);
    segment->vertex_map.Empty();
  }

  // Delete original faces
  for (int i = 0; i < original_faces.NEntries(); i++) {
    R3MeshFace *original_face = original_faces.Kth(i);
    mesh->DeleteFace(original_face);
  }

  // Delete original vertices
  for (int i = 0; i < original_vertices.NEntries(); i++) {
    R3MeshVertex *original_vertex = original_vertices.Kth(i);
    mesh->DeleteVertex(original_vertex);
  }

  // Process segments
  for (int i = 0; i < segments.NEntries(); i++) {
    Segment *segment = segments.Kth(i);

    // Gather info for vertices of segment
    RNArray<R3Point *> points;
    R3Vector sum_of_normals = R3zero_vector;
    for (int i = 0; i < segment->vertices.NEntries(); i++) {
      R3MeshVertex *vertex = segment->vertices.Kth(i);
      sum_of_normals += mesh->VertexNormal(vertex);
      points.Insert((R3Point *) &(mesh->VertexPosition(vertex)));
    }

    // Create a coordinate system for segment
    R3Point centroid = R3Centroid(points);
    R3Triad triad = R3PrincipleAxes(centroid, points);
    if (triad[2].Dot(sum_of_normals) < 0) triad.Rotate(triad[0], RN_PI);
    R3CoordSystem cs(centroid, triad);
    segment->matrix = cs.InverseMatrix();

    // SOMEHOW THE PARAMETERIZATION FOR 1/2 OF SEGMENTS IS FLIPPED BACKWARDS

    // Determine extent of segment in texture coordinates
    segment->bbox = R2null_box;
    for (int j = 0; j < segment->vertices.NEntries(); j++) {
      R3MeshVertex *vertex = segment->vertices.Kth(j);
      R3Point position = mesh->VertexPosition(vertex);
      position = segment->matrix * position;
      R2Point texcoords(position.X(), position.Y());
      segment->bbox.Union(texcoords);
    }
  }

  // Sort segments by ylength
  segments.Sort(CompareSegments);
  
  // Return success
  return 1;
}
    


static int 
DeleteSegments(R3Mesh *mesh, const RNArray<Segment *>& segments)
{
  // Delete each segment
  for (int i = 0; i < segments.NEntries(); i++) {
    Segment *segment = segments.Kth(i);
    delete segment;
  }

  // Return success
  return 1;
}



////////////////////////////////////////////////////////////////////////
// Parameterization stuff
////////////////////////////////////////////////////////////////////////

static int 
ParameterizeSegments(R3Mesh *mesh, const RNArray<Segment *>& segments)
{
  // Initialize location of segment in texture coordinates
  R2Point origin(0,0);
  RNLength max_xlength = 0;
  RNLength max_yorigin = mesh->BBox().LongestAxisLength();

  // Parameterize each segment
  for (int i = 0; i < segments.NEntries(); i++) {
    Segment *segment = segments.Kth(i);

    // Update matrix to start texcoords bbox for segment at origin
    R4Matrix matrix = segment->matrix;
    matrix[0][3] += origin.X() - segment->bbox.XMin();    
    matrix[1][3] += origin.Y() - segment->bbox.YMin();    
      
    // Assign texture coordinates to vertices
    for (int j = 0; j < segment->vertices.NEntries(); j++) {
      R3MeshVertex *vertex = segment->vertices.Kth(j);
      R3Point position = mesh->VertexPosition(vertex);
      position = matrix * position;
      R2Point texcoords(position.X(), position.Y());
      mesh->SetVertexTextureCoords(vertex, texcoords);
    }

    // Update max xlength
    if (segment->bbox.XLength() > max_xlength) max_xlength = segment->bbox.XLength();

    // Update origin
    origin[1] += segment->bbox.YLength();
    if (origin[1] > max_yorigin) {
      origin[0] += max_xlength;
      origin[1] = 0;
      max_xlength = 0;
    }
  }
  
  // Return success
  return 1;
}



static int
ParameterizeMesh(R3Mesh *mesh)
{
  // Start statistics
  RNTime start_time;
  start_time.Read();

  // Create segments
  RNArray<Segment *> segments;
  if (!CreateSegments(mesh, segments)) return 0;

  // Parameterize segments
  if (!ParameterizeSegments(mesh, segments)) return 0;
  
  // Delete segments
  if (!DeleteSegments(mesh, segments)) return 0;

  // Print statistics
  if (print_verbose) {
    printf("Parameterized mesh ...\n");
    printf("  Time = %.2f seconds\n", start_time.Elapsed());
    printf("  # Segments = %d\n", segments.NEntries());
    fflush(stdout);
  }

  // Return success
  return 1;
}




static int
FlattenMesh(R3Mesh *mesh)
{
  // This is useful for debugging
  if (!flatten) return 1;

  // Replace all vertex positions with texture coordinates
  for (int i = 0; i < mesh->NVertices(); i++) {
    R3MeshVertex *vertex = mesh->Vertex(i);
    R2Point texcoords = mesh->VertexTextureCoords(vertex);
    R3Point flattened_position(texcoords.X(), texcoords.Y(), 0);
    mesh->SetVertexPosition(vertex, flattened_position);
  }

  // Return success
  return 1;
}




////////////////////////////////////////////////////////////////////////
// Texture Creation Stuff
////////////////////////////////////////////////////////////////////////

static int
WriteTextures(R3Mesh *mesh, const char *directory_name, RNLength texel_spacing)
{
  // Check directory name
  if (!directory_name) return 1;
  
  // Compute bounding box of texture coordinates
  R2Box texcoords_bbox = R2null_box;
  for (int i = 0; i < mesh->NVertices(); i++) {
    R3MeshVertex *vertex = mesh->Vertex(i);
    const R2Point& texcoords = mesh->VertexTextureCoords(vertex);
    texcoords_bbox.Union(texcoords);
  }

  // Get useful variables for texcoords scaling
  R2Point t0, t1, t2;
  RNLength texcoords_xmin = texcoords_bbox.XMin();
  RNLength texcoords_ymin = texcoords_bbox.YMin();
  RNLength texcoords_xlength = texcoords_bbox.XLength();
  if (RNIsZero(texcoords_xlength)) return 0;
  RNLength texcoords_ylength = texcoords_bbox.YLength();
  if (RNIsZero(texcoords_ylength)) return 0;

  // Compute width and height from texel_spacing
  if (texel_spacing <= 0) texel_spacing = 1;
  int width = (int) (texcoords_bbox.XLength() / texel_spacing + 0.5);
  int height = (int) (texcoords_bbox.YLength() / texel_spacing + 0.5);
  if ((width <= 0) || (height == 0)) return 0;

  // Update bounding box from width and height 
  texcoords_bbox[1][0] = texcoords_xmin + texel_spacing*width;
  texcoords_bbox[1][1] = texcoords_ymin + texel_spacing*height;

  // Allocate images for texture channels (initialized with zeroes)
  R2Grid mask_channel(width, height, texcoords_bbox);
  R2Grid face_channel(mask_channel);
  R2Grid px_channel(mask_channel);
  R2Grid py_channel(mask_channel);
  R2Grid pz_channel(mask_channel);
  R2Grid nx_channel(mask_channel);
  R2Grid ny_channel(mask_channel);
  R2Grid nz_channel(mask_channel);
  R2Grid red_channel(mask_channel);
  R2Grid green_channel(mask_channel);
  R2Grid blue_channel(mask_channel);
  R2Grid curvature_channel(mask_channel);
  R2Grid category_channel(mask_channel);
  R2Grid segment_channel(mask_channel);
  R2Grid material_channel(mask_channel);

  // Rasterize mesh face properties into texture channels
  for (int i = 0; i < mesh->NFaces(); i++) {
    R3MeshFace *face = mesh->Face(i);
    int category = mesh->FaceCategory(face);
    int segment = mesh->FaceSegment(face);
    int material = mesh->FaceMaterial(face);
    R3Vector normal = mesh->FaceNormal(face);
    R3MeshVertex *v0 = mesh->VertexOnFace(face, 0);
    R3MeshVertex *v1 = mesh->VertexOnFace(face, 1);
    R3MeshVertex *v2 = mesh->VertexOnFace(face, 2);
    RNRgb rgb0 = mesh->VertexColor(v0);
    RNRgb rgb1 = mesh->VertexColor(v1);
    RNRgb rgb2 = mesh->VertexColor(v2);
    RNScalar curvature0 = mesh->VertexMeanCurvature(v0);
    RNScalar curvature1 = mesh->VertexMeanCurvature(v1);
    RNScalar curvature2 = mesh->VertexMeanCurvature(v2);
    R3Point p0 = mesh->VertexPosition(v0);
    R3Point p1 = mesh->VertexPosition(v1);
    R3Point p2 = mesh->VertexPosition(v2);
    R2Point uv0 = mesh->VertexTextureCoords(v0);
    R2Point uv1 = mesh->VertexTextureCoords(v1);
    R2Point uv2 = mesh->VertexTextureCoords(v2);
    t0[0] = width * (uv0[0] - texcoords_xmin) / texcoords_xlength;
    t0[1] = height * (uv0[1] - texcoords_ymin) / texcoords_ylength;
    t1[0] = width * (uv1[0] - texcoords_xmin) / texcoords_xlength;
    t1[1] = height * (uv1[1] - texcoords_ymin) / texcoords_ylength;
    t2[0] = width * (uv2[0] - texcoords_xmin) / texcoords_xlength;
    t2[1] = height * (uv2[1] - texcoords_ymin) / texcoords_ylength;
    mask_channel.RasterizeGridTriangle(t0, t1, t2, 1.0, R2_GRID_REPLACE_OPERATION);
    face_channel.RasterizeGridTriangle(t0, t1, t2, i+1, R2_GRID_REPLACE_OPERATION);
    category_channel.RasterizeGridTriangle(t0, t1, t2, category, category, category, R2_GRID_REPLACE_OPERATION);
    segment_channel.RasterizeGridTriangle(t0, t1, t2, segment, segment, segment, R2_GRID_REPLACE_OPERATION);
    material_channel.RasterizeGridTriangle(t0, t1, t2, material, material, material, R2_GRID_REPLACE_OPERATION);
    nx_channel.RasterizeGridTriangle(t0, t1, t2, normal.X(), R2_GRID_REPLACE_OPERATION);
    ny_channel.RasterizeGridTriangle(t0, t1, t2, normal.Y(), R2_GRID_REPLACE_OPERATION);
    nz_channel.RasterizeGridTriangle(t0, t1, t2, normal.Z(), R2_GRID_REPLACE_OPERATION);
    px_channel.RasterizeGridTriangle(t0, t1, t2, p0.X(), p1.X(), p2.X(), R2_GRID_REPLACE_OPERATION);
    py_channel.RasterizeGridTriangle(t0, t1, t2, p0.Y(), p1.Y(), p2.Y(), R2_GRID_REPLACE_OPERATION);
    pz_channel.RasterizeGridTriangle(t0, t1, t2, p0.Z(), p1.Z(), p2.Z(), R2_GRID_REPLACE_OPERATION);
    red_channel.RasterizeGridTriangle(t0, t1, t2, rgb0.R(), rgb1.R(), rgb2.R(), R2_GRID_REPLACE_OPERATION);
    green_channel.RasterizeGridTriangle(t0, t1, t2, rgb0.G(), rgb1.G(), rgb2.G(), R2_GRID_REPLACE_OPERATION);
    blue_channel.RasterizeGridTriangle(t0, t1, t2, rgb0.B(), rgb1.B(), rgb2.B(), R2_GRID_REPLACE_OPERATION);
    curvature_channel.RasterizeGridTriangle(t0, t1, t2, curvature0, curvature1, curvature2, R2_GRID_REPLACE_OPERATION);
  }

  // Create color image (combined 3 channels into a single png)
  RNRgb color(0,0,0);
  R2Image color_image(width, height, 3);
  for (int i = 0; i < width; i++) {
    for (int j = 0; j < height; j++) {
      color[0] = red_channel.GridValue(i, j);
      color[1] = green_channel.GridValue(i, j);
      color[2] = blue_channel.GridValue(i, j);
      color_image.SetPixelRGB(i, j, color);
    }
  }

  // Clamp category, segment, and material channels (for png output)
  category_channel.Threshold(0, 0, R2_GRID_KEEP_VALUE);
  segment_channel.Threshold(0, 0, R2_GRID_KEEP_VALUE);
  material_channel.Threshold(0, 0, R2_GRID_KEEP_VALUE);
  category_channel.Threshold(65536, R2_GRID_KEEP_VALUE, 65536);
  segment_channel.Threshold(65536, R2_GRID_KEEP_VALUE, 65536);
  material_channel.Threshold(65536, R2_GRID_KEEP_VALUE, 65536);

  // Scale normal channels = 65535 * (value+1.0)/2.0 (for png output)
  nx_channel.Add(1.0); nx_channel.Multiply(65535/2.0);
  ny_channel.Add(1.0); ny_channel.Multiply(65535/2.0);
  nz_channel.Add(1.0); nz_channel.Multiply(65535/2.0);
  nx_channel.Threshold(0, 0, R2_GRID_KEEP_VALUE);
  ny_channel.Threshold(0, 0, R2_GRID_KEEP_VALUE);
  nz_channel.Threshold(0, 0, R2_GRID_KEEP_VALUE);
  nx_channel.Threshold(65535, R2_GRID_KEEP_VALUE, 65535);
  ny_channel.Threshold(65535, R2_GRID_KEEP_VALUE, 65535);
  nz_channel.Threshold(65535, R2_GRID_KEEP_VALUE, 65535);

  // Scale position channel (for png output)
  px_channel.Multiply(1000);
  py_channel.Multiply(1000);
  pz_channel.Multiply(1000);
  px_channel.Threshold(0, 0, R2_GRID_KEEP_VALUE);
  py_channel.Threshold(0, 0, R2_GRID_KEEP_VALUE);
  pz_channel.Threshold(0, 0, R2_GRID_KEEP_VALUE);
  px_channel.Threshold(65535, R2_GRID_KEEP_VALUE, 65535);
  py_channel.Threshold(65535, R2_GRID_KEEP_VALUE, 65535);
  pz_channel.Threshold(65535, R2_GRID_KEEP_VALUE, 65535);

  // Scale curvature channel = 65535 * (value+100)/2.0 (for png output)
  curvature_channel.Threshold(-100, -100, R2_GRID_KEEP_VALUE);
  curvature_channel.Threshold(100, R2_GRID_KEEP_VALUE, 100);
  curvature_channel.Add(100);
  curvature_channel.Multiply(65535/200.0);
  curvature_channel.Threshold(0, 0, R2_GRID_KEEP_VALUE);
  curvature_channel.Threshold(65535, R2_GRID_KEEP_VALUE, 65535);

  // Create texture directory
  char cmd[1024];
  sprintf(cmd, "mkdir -p %s", directory_name);
  system(cmd);

  // Write channels
  char filename[1024];
  sprintf(filename, "%s/mask.png", directory_name);
  if (!mask_channel.WriteFile(filename)) return 0;
  sprintf(filename, "%s/face.pfm", directory_name);
  if (!face_channel.WriteFile(filename)) return 0;
  sprintf(filename, "%s/px.png", directory_name);
  if (!px_channel.WriteFile(filename)) return 0;
  sprintf(filename, "%s/py.png", directory_name);
  if (!py_channel.WriteFile(filename)) return 0;
  sprintf(filename, "%s/pz.png", directory_name);
  if (!pz_channel.WriteFile(filename)) return 0;
  sprintf(filename, "%s/nx.png", directory_name);
  if (!nx_channel.WriteFile(filename)) return 0;
  sprintf(filename, "%s/ny.png", directory_name);
  if (!ny_channel.WriteFile(filename)) return 0;
  sprintf(filename, "%s/nz.png", directory_name);
  if (!nz_channel.WriteFile(filename)) return 0;
  sprintf(filename, "%s/curvature.png", directory_name);
  if (!curvature_channel.WriteFile(filename)) return 0;
  sprintf(filename, "%s/category.png", directory_name);
  if (!category_channel.WriteFile(filename)) return 0;
  sprintf(filename, "%s/segment.png", directory_name);
  if (!segment_channel.WriteFile(filename)) return 0;
  sprintf(filename, "%s/material.png", directory_name);
  if (!material_channel.WriteFile(filename)) return 0;
  sprintf(filename, "%s/color.png", directory_name);
  if (!color_image.Write(filename)) return 0;

  // Return success
  return 1;
}



////////////////////////////////////////////////////////////////////////
// Argument Parsing Stuff
////////////////////////////////////////////////////////////////////////

static int 
ParseArgs(int argc, char **argv)
{
  // Parse arguments
  argc--; argv++;
  while (argc > 0) {
    if ((*argv)[0] == '-') {
      if (!strcmp(*argv, "-v")) print_verbose = 1;
      else if (!strcmp(*argv, "-flatten")) flatten = 1;
      else if (!strcmp(*argv, "-separate_face_segments")) separate_face_segments = 1;
      else if (!strcmp(*argv, "-separate_face_materials")) separate_face_materials = 1;
      else if (!strcmp(*argv, "-output_textures")) { argc--; argv++; output_texture_directory = *argv; }
      else if (!strcmp(*argv, "-texel_spacing")) { argc--; argv++; texel_spacing = atof(*argv); }
      else { fprintf(stderr, "Invalid program argument: %s", *argv); exit(1); }
      argv++; argc--;
    }
    else {
      if (!input_mesh_filename) input_mesh_filename = *argv;
      else if (!output_mesh_filename) output_mesh_filename = *argv;
      else { fprintf(stderr, "Invalid program argument: %s", *argv); exit(1); }
      argv++; argc--;
    }
  }

  // Check filenames
  if (!input_mesh_filename || !output_mesh_filename) {
    fprintf(stderr, "Usage: mshparam inputmesh outputmesh\n");
    return 0;
  }

  // Return OK status 
  return 1;
}



////////////////////////////////////////////////////////////////////////
// Main
////////////////////////////////////////////////////////////////////////

int 
main(int argc, char **argv)
{
  // Check number of arguments
  if (!ParseArgs(argc, argv)) exit(1);

  // Read input mesh
  R3Mesh *mesh = ReadMesh(input_mesh_filename);
  if (!mesh) exit(-1);

  // Segment mesh
  if (!ParameterizeMesh(mesh)) exit(-1);

  // Flatten mesh
  if (!FlattenMesh(mesh)) exit(-1);

  // Write mesh
  if (!WriteMesh(mesh, output_mesh_filename)) exit(-1);

  // Write textures
  if (!WriteTextures(mesh, output_texture_directory, texel_spacing)) exit(-1);

  // Return success 
  return 0;
}



