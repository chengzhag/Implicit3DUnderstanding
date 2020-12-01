// Program to generate points set from mesh



////////////////////////////////////////////////////////////////////////
// Include files
////////////////////////////////////////////////////////////////////////

namespace gaps {}
using namespace gaps;
#include "R3Shapes/R3Shapes.h"



////////////////////////////////////////////////////////////////////////
// Constant definitions
////////////////////////////////////////////////////////////////////////

enum {
  RANDOM_VERTICES,
  RANDOM_SURFACE_POINTS,
  ITERATIVE_FURTHEST_VERTEX,  
  EXTERIOR_SURFACE_POINTS,
  NEAR_SURFACE_POINTS,
  VISIBLE_SURFACE_POINTS,
  UNIFORM_IN_BBOX,
  SCALE_SPACE_MAXIMA,
  SCALE_SPACE_MINIMA,
  SCALE_SPACE_EXTREMA,
  PROPERTY_MAXIMA,
  PROPERTY_MINIMA,
  PROPERTY_EXTREMA,
  CENTER_OF_MASS,
  VERTEX_CLOSEST_TO_CENTER_OF_MASS,
  NUM_SELECTION_METHODS
};



////////////////////////////////////////////////////////////////////////
// Program arguments
////////////////////////////////////////////////////////////////////////

static const char *mesh_name  = NULL;
static const char *points_name = NULL;
static const char *property_name  = NULL;
static int selection_method = RANDOM_SURFACE_POINTS;
static int min_points = -1;
static int max_points = -1;
static int num_points = -1;
static R3Box bbox(FLT_MAX, FLT_MAX, FLT_MAX, -FLT_MAX, -FLT_MAX, -FLT_MAX);
static double min_spacing = -1;
static double min_normalized_spacing = -1;
static double min_relative_spacing = -1;
static double max_distance = FLT_MAX;
static double near_surface_bias_exponent = 2;
static double curvature_exponent = 0;
static RNScalar curvature_max = 100;
static RNBoolean binary_sdf = FALSE;
static RNBoolean print_verbose = FALSE;
static RNBoolean print_debug = FALSE;



////////////////////////////////////////////////////////////////////////
// Type definitions
////////////////////////////////////////////////////////////////////////

struct Point {
  Point(void) 
    : vertex(NULL), face(NULL), position(0,0,0), normal(0,0,0), sdf(0) {};
  Point(const R3Point& position, const R3Vector& normal) 
    : vertex(NULL), face(NULL), position(position), normal(normal), sdf(0) {};
  Point(R3Mesh *mesh, R3MeshVertex *vertex)
    : vertex(vertex), face(NULL), position(mesh->VertexPosition(vertex)), normal(mesh->VertexNormal(vertex)), sdf(0) {};
  R3MeshVertex *vertex;
  R3MeshFace *face;
  R3Point position;
  R3Vector normal;
  RNScalar sdf;
};



////////////////////////////////////////////////////////////////////////
// Mesh file input
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



static R3MeshProperty *
ReadProperty(R3Mesh *mesh, const char *filename)
{
  // Start statistics
  RNTime start_time;
  start_time.Read();

  // Allocate properties
  R3MeshProperty *property = new R3MeshProperty(mesh);
  if (!property) {
    RNFail("Unable to allocate property for %s\n", filename);
    return NULL;
  }

  // Read property from file
  if (!property->Read(filename)) {
    delete property;
    return NULL;
  }

  // Print statistics
  if (print_verbose) {
    printf("Read property from %s ...\n", filename);
    printf("  Time = %.2f seconds\n", start_time.Elapsed());
    printf("  Minimum = %g\n", property->Minimum());
    printf("  Maximum = %g\n", property->Maximum());
    printf("  Median = %g\n", property->Median());
    printf("  Mean = %g\n", property->Mean());
    printf("  # Minima = %d\n", property->LocalMinimumCount());
    printf("  # Maxima = %d\n", property->LocalMaximumCount());
    printf("  L1Norm = %g\n", property->L1Norm());
    printf("  L2Norm = %g\n", property->L2Norm());
    fflush(stdout);
  }

  // Return property 
  return property;
}



////////////////////////////////////////////////////////////////////////
// Point file output
////////////////////////////////////////////////////////////////////////

RNBoolean 
WritePoints(R3Mesh *mesh, const RNArray<Point *>& points, const char *filename)
{
  // Start statistics
  RNTime start_time;
  start_time.Read();

  // Parse input filename extension
  const char *extension;
  if (!(extension = strrchr(filename, '.'))) {
    printf("Filename %s has no extension (e.g., .ply)\n", filename);
    return 0;
  }

  // Check file type
  if (!strcmp(extension, ".pid")) {
    // Open file
    FILE *fp = stdout;
    if (filename) {
      fp = fopen(filename, "w");
      if (!fp) {
        RNFail("Unable to open output file %s\n", filename);
        return FALSE;
      }
    }

    // Write vertex ids
    for (int i = 0; i < points.NEntries(); i++) {
      Point *point = points[i];
      if (!point->vertex) { RNFail("Incompatible output type: pid\n"); return 0; }
      int vertex_id = (point->vertex) ? mesh->VertexID(point->vertex) : -1;
      fprintf(fp, "%d\n", vertex_id);
    }

    // Close file
    if (filename) fclose(fp);
  }
  else if (!strcmp(extension, ".vts")) {
    // Open file
    FILE *fp = stdout;
    if (filename) {
      fp = fopen(filename, "w");
      if (!fp) {
        RNFail("Unable to open output file %s\n", filename);
        return FALSE;
      }
    }

    // Write points
    for (int i = 0; i < points.NEntries(); i++) {
      Point *point = points[i];
      int vertex_id = (point->vertex) ? mesh->VertexID(point->vertex) : -1;
      const R3Point& position = point->position;
      fprintf(fp, "%d %g %g %g\n", vertex_id, position.X(), position.Y(), position.Z());
    }

    // Close file
    if (filename) fclose(fp);
  }
  else if (!strcmp(extension, ".xyzn")) {
    // Open file
    FILE *fp = stdout;
    if (filename) {
      fp = fopen(filename, "w");
      if (!fp) {
        RNFail("Unable to open output file %s\n", filename);
        return FALSE;
      }
    }

    // Write points
    for (int i = 0; i < points.NEntries(); i++) {
      Point *point = points[i];
      fprintf(fp, "%g %g %g   %g %g %g\n", 
        point->position.X(), point->position.Y(), point->position.Z(),
        point->normal.X(), point->normal.Y(), point->normal.Z());
    }

    // Close file
    if (filename) fclose(fp);
  }
  else if (!strcmp(extension, ".xyz")) {
    // Open file
    FILE *fp = stdout;
    if (filename) {
      fp = fopen(filename, "w");
      if (!fp) {
        RNFail("Unable to open output file %s\n", filename);
        return FALSE;
      }
    }

    // Write points
    for (int i = 0; i < points.NEntries(); i++) {
      Point *point = points[i];
      fprintf(fp, "%g %g %g\n", point->position.X(), point->position.Y(), point->position.Z());
    }

    // Close file
    if (filename) fclose(fp);
  }
  else if (!strcmp(extension, ".pts")) {
    // Open file
    FILE *fp = stdout;
    if (filename) {
      fp = fopen(filename, "wb");
      if (!fp) {
        RNFail("Unable to open output file %s\n", filename);
        return FALSE;
      }
    }

    // Write points
    float coordinates[6];
    for (int i = 0; i < points.NEntries(); i++) {
      Point *point = points[i];
      coordinates[0] = point->position.X();
      coordinates[1] = point->position.Y();
      coordinates[2] = point->position.Z();
      coordinates[3] = point->normal.X();
      coordinates[4] = point->normal.Y();
      coordinates[5] = point->normal.Z();
      if (fwrite(coordinates, sizeof(float), 6, fp) != (unsigned int) 6) {
        RNFail("Unable to write point to output file %s\n", filename);
        return FALSE;
      }
    }

    // Close file
    if (filename) fclose(fp);
  }
  else if (!strcmp(extension, ".sdf")) {
    // Open binary file
    FILE *fp = fopen(filename, "wb");
    if (!fp) {
      RNFail("Unable to open output file %s\n", filename);
      return 0;
    }

    // Write points
    float values[4];
    for (int i = 0; i < points.NEntries(); i++) {
      Point *point = points.Kth(i);
      values[0] = point->position.X();
      values[1] = point->position.Y();
      values[2] = point->position.Z();
      values[3] = point->sdf;
      if (fwrite(values, sizeof(float), 4, fp) != (unsigned int) 4) {
        RNFail("Unable to write point to output file %s\n", filename);
        return 0;
      }
    }
  }
  else {
    // Create mesh
    R3Mesh mesh;
    for (int i = 0; i < points.NEntries(); i++) {
      Point *point = points[i];
      mesh.CreateVertex(point->position, point->normal);
    }

    // Write mesh file
    if (!mesh.WriteFile(filename)) return FALSE;
  }

  // Print message
  if (print_verbose) {
    fprintf(stdout, "Wrote points to %s ...\n", (filename) ? filename : "stdout");
    printf("  Time = %.2f seconds\n", start_time.Elapsed());
    fflush(stdout);
  }

  // Return success
  return TRUE;
}



////////////////////////////////////////////////////////////////////////
// Sample weighting utility functions
////////////////////////////////////////////////////////////////////////

static RNScalar
ComputeFaceSamplingWeights(R3Mesh *mesh, R3MeshProperty *property,
  RNScalar curvature_exponent = 0, RNScalar curvature_max = RN_INFINITY)
{
  // Compute face weights
  RNScalar total_weight = 0;
  for (int i = 0; i < mesh->NFaces(); i++) {
    R3MeshFace *face = mesh->Face(i);

    // Initailize weight by face area
    RNScalar weight = mesh->FaceArea(face);

    // Scale weight by curvatures
    if (property) {
      // Compute average value of vertices
      R3MeshVertex *v0 = mesh->VertexOnFace(face, 0);
      R3MeshVertex *v1 = mesh->VertexOnFace(face, 1);
      R3MeshVertex *v2 = mesh->VertexOnFace(face, 2);
      RNScalar c0 = fabs(property->VertexValue(v0));
      RNScalar c1 = fabs(property->VertexValue(v1));
      RNScalar c2 = fabs(property->VertexValue(v2));
      weight *= 1.0 + (c0 + c1 + c2) / 3.0;
    }
    else if (curvature_exponent > 0) {
      // Compute curvature
      R3MeshVertex *v0 = mesh->VertexOnFace(face, 0);
      R3MeshVertex *v1 = mesh->VertexOnFace(face, 1);
      R3MeshVertex *v2 = mesh->VertexOnFace(face, 2);
      RNScalar c0 = fabs(mesh->VertexMeanCurvature(v0));
      RNScalar c1 = fabs(mesh->VertexMeanCurvature(v1));
      RNScalar c2 = fabs(mesh->VertexMeanCurvature(v2));
      RNScalar curvature = 0;
      curvature += (c0 < curvature_max) ? c0 : curvature_max;
      curvature += (c1 < curvature_max) ? c1 : curvature_max;
      curvature += (c2 < curvature_max) ? c2 : curvature_max;
      curvature /= 3.0;

      // Adjust weight
      weight *= pow(1.0 + curvature, curvature_exponent);
    }

    // Set face sampling weight
    mesh->SetFaceValue(face, weight);
    total_weight += weight;
  }

  // Delete property
  if (property) delete property;

  // Return total weight
  return total_weight;
}



////////////////////////////////////////////////////////////////////////
// Signed distance utility functions
////////////////////////////////////////////////////////////////////////

static R3Grid *
CreateSignedDistanceGrid(R3Mesh *mesh, const R3Box& b)
{
  // Allocate grid
  int grid_resolution = 256;
  RNScalar target_spacing = b.LongestAxisLength() / grid_resolution;

  // Allocate grid
  R3Grid *grid = new R3Grid(b, target_spacing, 5, 2*grid_resolution, 5);
  if (!grid) {
    RNFail("Unable to allocate grid\n");
    return NULL;
  }

  // Rasterize each triangle into grid
  for (int i = 0; i < mesh->NFaces(); i++) {
    R3MeshFace *face = mesh->Face(i);
    const R3Point& p0 = mesh->VertexPosition(mesh->VertexOnFace(face, 0));
    const R3Point& p1 = mesh->VertexPosition(mesh->VertexOnFace(face, 1));
    const R3Point& p2 = mesh->VertexPosition(mesh->VertexOnFace(face, 2));
    grid->RasterizeWorldTriangle(p0, p1, p2, i+1, R3_GRID_REPLACE_OPERATION);
  }

  // Estimate distance to closest face (at grid resolution)
  grid->SquaredDistanceTransform();
  grid->Sqrt();
  grid->Multiply(grid->GridToWorldScaleFactor());
  
  // Initialize sign grid
  R3Grid sign_grid(*grid);
  sign_grid.Clear(-1);

  // Seed search with border voxels (assumed outside)
  RNArray<const RNScalar *> stack;
  const RNScalar *sign_grid_values = sign_grid.GridValues();
  RNLength grid_spacing = grid->GridToWorldScaleFactor();
  RNScalar max_grid_value = sqrt(3)*grid_spacing;
  int ix, iy, iz, grid_index;
  for (int i = 0; i < sign_grid.XResolution(); i++) {
    for (int j = 0; j < sign_grid.YResolution(); j++) {
      for (int k = 0; k < sign_grid.ZResolution(); k++) {
        if ((i == 0) || (i == sign_grid.XResolution()-1) ||
            (j == 0) || (j == sign_grid.YResolution()-1) ||
            (k == 0) || (k == sign_grid.ZResolution()-1)) {
          if (grid->GridValue(i, j, k) < grid_spacing) continue;
          sign_grid.IndicesToIndex(i, j, k, grid_index);
          stack.Insert(&sign_grid_values[grid_index]);
          sign_grid.SetGridValue(grid_index, 1);
        }
      }
    }
  }

  // Flood fill
  while (!stack.IsEmpty()) {
    const RNScalar *sign_grid_valuep = stack.Tail(); stack.RemoveTail();
    int grid_index = sign_grid_valuep - sign_grid_values;
    assert((grid_index >= 0) && (grid_index < sign_grid.NEntries()));
    sign_grid.IndexToIndices(grid_index, ix, iy, iz);
    for (int dz = -1; dz <= 1; dz++) {
      int gz = iz + dz;
      if ((gz < 0) || (gz >= sign_grid.ZResolution())) continue;
      for (int dy = -1; dy <= 1; dy++) {
        int gy = iy + dy;
        if ((gy < 0) || (gy >= sign_grid.YResolution())) continue;
        for (int dx = -1; dx <= 1; dx++) {
          int gx = ix + dx;
          if ((gx < 0) || (gx >= sign_grid.XResolution())) continue;
          if (sign_grid.GridValue(gx, gy, gz) > 0) continue;

          // Check if reached the boundary
          RNScalar d = grid->GridValue(ix, iy, iz);
          if (d < max_grid_value) continue;
             
          // Mark cell as outside and continue DFS 
          sign_grid.IndicesToIndex(gx, gy, gz, grid_index);
          stack.Insert(&sign_grid_values[grid_index]);
          sign_grid.SetGridValue(grid_index, 1);
        }
      }
    }
  }
  
  // Apply signs
  sign_grid.Threshold(0, -1, 1);
  grid->Multiply(sign_grid);

  // Return signed distance grid
  return grid;
}



////////////////////////////////////////////////////////////////////////
// Random surface point sampling
////////////////////////////////////////////////////////////////////////

static RNArray<Point *> *
SelectRandomSurfacePoints(R3Mesh *mesh, R3MeshProperty *property, int npoints)
{
  // Allocate array of points
  RNArray<Point *> *points = new RNArray<Point *>();
  if (!points) {
    RNFail("Unable to allocate array of points\n");
    return NULL;
  }

  // Count total value of faces
  RNArea total_value = ComputeFaceSamplingWeights(mesh,
    property, curvature_exponent, curvature_max);
  if (total_value == 0) return points;
  
  // Generate points
  RNSeedRandomScalar();
  for (int i = 0; i < mesh->NFaces(); i++) {
    R3MeshFace *face = mesh->Face(i);

    // Get vertex positions
    R3MeshVertex *v0 = mesh->VertexOnFace(face, 0);
    R3MeshVertex *v1 = mesh->VertexOnFace(face, 1);
    R3MeshVertex *v2 = mesh->VertexOnFace(face, 2);
    const R3Point& p0 = mesh->VertexPosition(v0);
    const R3Point& p1 = mesh->VertexPosition(v1);
    const R3Point& p2 = mesh->VertexPosition(v2);
    const R3Vector& n0 = mesh->VertexNormal(v0);
    const R3Vector& n1 = mesh->VertexNormal(v1);
    const R3Vector& n2 = mesh->VertexNormal(v2);

    // Determine number of points for face 
    RNScalar ideal_face_npoints = npoints * mesh->FaceValue(face) / total_value;
    int face_npoints = (int) ideal_face_npoints;
    RNScalar remainder = ideal_face_npoints - face_npoints;
    if (remainder > RNRandomScalar()) face_npoints++;

    // Generate random points in face
    for (int j = 0; j < face_npoints; j++) {
      RNScalar r1 = sqrt(RNRandomScalar());
      RNScalar r2 = RNRandomScalar();
      RNScalar t0 = (1.0 - r1);
      RNScalar t1 = r1 * (1.0 - r2);
      RNScalar t2 = r1 * r2;
      R3Point position = t0*p0 + t1*p1 + t2*p2;
      R3Vector normal = t0*n0 + t1*n1 + t2*n2; normal.Normalize();
      Point *point = new Point(position, normal);
      points->Insert(point);
    }
  }

  // Return points
  return points;
}



////////////////////////////////////////////////////////////////////////
// Vertex sampling
////////////////////////////////////////////////////////////////////////

static RNArray<Point *> *
SelectVertexPoints(R3Mesh *mesh, int npoints, double min_spacing)
{
  // Check number of points
  if (npoints > mesh->NVertices()) npoints = mesh->NVertices();

  // Allocate array of points
  RNArray<Point *> *points = new RNArray<Point *>();
  if (!points) {
    RNFail("Unable to allocate array of points\n");
    return NULL;
  }

  // Set vertex values = surface area associated with vertex
  RNScalar total_value = 0;
  for (int i = 0; i < mesh->NVertices(); i++) {
    R3MeshVertex *vertex = mesh->Vertex(i);
    RNScalar value = mesh->VertexArea(vertex);
    mesh->SetVertexValue(vertex, value);
    total_value += value;
  }

  // Select vertices weighted by value
  R3mesh_mark++;
  int max_iterations = 100 * npoints;
  for (int i = 0; (npoints > 0) && (i < max_iterations); i++) {
    RNScalar cumulative_value = 0;
    RNScalar target_value = RNRandomScalar() * total_value;
    for (int j = 0; j < mesh->NVertices(); j++) {
      R3MeshVertex *vertex = mesh->Vertex(j);
      if (mesh->VertexMark(vertex) != R3mesh_mark) {
        RNScalar vertex_value = mesh->VertexValue(vertex);
        cumulative_value += vertex_value;
        if ((cumulative_value >= target_value) || (j == mesh->NVertices()-1)) {
          // Insert point
          Point *point = new Point(mesh, vertex);
          points->Insert(point);
          npoints--;

          // Mark point 
          mesh->SetVertexMark(vertex, R3mesh_mark);
          total_value -= vertex_value;

          // Mark nearby points
          if (min_spacing > 0) {
            RNLength *distances = mesh->DijkstraDistances(vertex, min_spacing);
            for (int k = 0; k < mesh->NVertices(); k++) {
              R3MeshVertex *nearby_vertex = mesh->Vertex(k);
              if (mesh->VertexMark(nearby_vertex) != R3mesh_mark) {
                if (distances[k] <= min_spacing) {
                  // Mark point within min_spacing
                  mesh->SetVertexMark(nearby_vertex, R3mesh_mark);
                  total_value -= vertex_value;
                }
              }
            }
            delete [] distances;
          }

          // Break
          break;
        }
      }
    }
  }

  // Return points
  return points;
}



////////////////////////////////////////////////////////////////////////
// Furthest surface point sampling
////////////////////////////////////////////////////////////////////////

// Type definitions

struct VertexData {
  R3MeshVertex *vertex;
  double distance;
  R3MeshVertex *closest_seed;
  VertexData **heappointer;
};



static R3MeshVertex *
FindFurthestVertex(R3Mesh *mesh, const RNArray<R3MeshVertex *>& seeds, VertexData *vertex_data)
{
  // Initialize vertex data
  for (int i = 0; i < mesh->NVertices(); i++) {
    VertexData *data = &vertex_data[i];
    data->distance = FLT_MAX;
    data->closest_seed = NULL;
    data->heappointer = NULL;
  }

  // Initialize heap
  VertexData tmp;
  RNHeap<VertexData *> heap(&tmp, &(tmp.distance), &(tmp.heappointer));
  for (int i = 0; i < seeds.NEntries(); i++) {
    R3MeshVertex *vertex = seeds.Kth(i);
    VertexData *data = &vertex_data[mesh->VertexID(vertex)];
    data->closest_seed = vertex;
    data->distance = 0;
    heap.Push(data);
  }

  // Iteratively pop off heap until find furthest vertex
  R3MeshVertex *vertex = NULL;
  while (!heap.IsEmpty()) {
    VertexData *data = heap.Pop();
    vertex = data->vertex;
    for (int j = 0; j < mesh->VertexValence(vertex); j++) {
      R3MeshEdge *edge = mesh->EdgeOnVertex(vertex, j);
      R3MeshVertex *neighbor_vertex = mesh->VertexAcrossEdge(edge, vertex);
      VertexData *neighbor_data = (VertexData *) mesh->VertexData(neighbor_vertex);
      RNScalar old_distance = neighbor_data->distance;
      RNScalar new_distance = mesh->EdgeLength(edge) + data->distance;
      if (new_distance < old_distance) {
        neighbor_data->distance = new_distance;
        neighbor_data->closest_seed = data->closest_seed;
        if (old_distance < FLT_MAX) heap.Update(neighbor_data);
        else heap.Push(neighbor_data);
      }
    }
  }

  // Return furthest vertex
  return vertex;
}



static RNArray<Point *> *
SelectFurthestPoints(R3Mesh *mesh, const RNArray<R3MeshVertex *>& seeds, int npoints, double min_spacing)
{
  // Check/adjust number of points
  if (npoints > mesh->NVertices()) npoints = mesh->NVertices();

  // Allocate array of points
  RNArray<Point *> *points = new RNArray<Point *>();
  if (!points) {
    RNFail("Unable to allocate array of points\n");
    return NULL;
  }

  // Allocate vertex data
  VertexData *vertex_data = new VertexData [ mesh->NVertices() ];
  if (!vertex_data) {
    RNFail("Unable to allocate vertex data\n");
    return NULL;
  }

  // Initialize vertex data
  for (int i = 0; i < mesh->NVertices(); i++) {
    R3MeshVertex *vertex = mesh->Vertex(i);
    VertexData *data = &vertex_data[i];
    mesh->SetVertexData(vertex, data);
    data->vertex = vertex;
    data->distance = FLT_MAX;
    data->closest_seed = NULL;
    data->heappointer = NULL;
  }

  // Copy seeds 
  RNArray<R3MeshVertex *> vertices;
  for (int i = 0; i < seeds.NEntries(); i++) {
    R3MeshVertex *vertex = seeds.Kth(i);
    Point *point = new Point(mesh, vertex);
    points->Insert(point);
    vertices.Insert(vertex);
  }

  // Create random starting vertex in each connected component, if none provided
  RNArray<R3MeshVertex *> tmp;
  if (seeds.IsEmpty()) {
    int *marks = new int [ mesh->NVertices() ];
    for (int i = 0; i < mesh->NVertices(); i++) marks[i] = 0;
    for (int i = 0; i < mesh->NVertices(); i++) {
      R3MeshVertex *vertex = mesh->Vertex(i);

      // Check if already marked
      if (marks[i] != 0) continue;
      marks[i] = 1;

      // Check if completely disconnected
      if (mesh->VertexValence(vertex) == 0) continue;

      // Select vertex as temporary seed
      vertices.Insert(vertex);
      tmp.Insert(vertex);

      // Mark connected component
      RNArray<R3MeshVertex *> connected_component;
      mesh->FindConnectedVertices(vertex, connected_component);
      for (int j = 0; j < connected_component.NEntries(); j++) {
        R3MeshVertex *connected_vertex = connected_component.Kth(j);
        marks[mesh->VertexID(connected_vertex)] = 1;
      }
    }
  }

  // Iteratively find furthest vertex
  while (vertices.NEntries() - tmp.NEntries() < npoints) {
    R3MeshVertex *vertex = FindFurthestVertex(mesh, vertices, vertex_data);
    if (!vertex) break;
    VertexData *data = &vertex_data[mesh->VertexID(vertex)];
    if ((min_spacing > 0) && (data->distance < min_spacing)) break;
    if (tmp.FindEntry(data->closest_seed)) { 
      tmp.Remove(data->closest_seed); 
      vertices.Remove(data->closest_seed);
    }
    Point *point = new Point(mesh, vertex);
    points->Insert(point);
    vertices.Insert(vertex);
  }

  // Delete vertex data
  delete [] vertex_data;

  // Return points
  return points;
}



#if 0
static RNArray<Point *> *
SelectFurthestPoints(R3Mesh *mesh, const RNArray<Point *>& seeds, int npoints, double min_spacing)
{
  // Seeed with vertices closest to points
  RNArray<R3MeshVertex *> seed_vertices;
  for (int i = 0; i < seeds.NEntries(); i++) {
    R3MeshVertex *vertex = seeds[i]->vertex;
    if (vertex) seed_vertices.Insert(vertex);
    else { RNFail("Not implemented\n"); return NULL; }
  }

  // Select furthest points
  return SelectFurthestPoints(mesh, seed_vertices, npoints, min_spacing);
} 
#endif



static RNArray<Point *> *
SelectFurthestPoints(R3Mesh *mesh, int npoints, double min_spacing)
{
  // Select furthest points from scratch
  RNArray<R3MeshVertex *> seeds;
  return SelectFurthestPoints(mesh, seeds, npoints, min_spacing);
}



////////////////////////////////////////////////////////////////////////
// Property extrema sampling
////////////////////////////////////////////////////////////////////////

static RNBoolean 
IsLocalExtremum(R3Mesh *mesh, R3MeshVertex *vertex, R3MeshProperty *property, int extrema_type, RNScalar epsilon = 0)
{
  // Return whether value at vertex is local extremum
  RNBoolean minimum = TRUE;
  RNBoolean maximum = TRUE;
  RNScalar value = property->VertexValue(mesh->VertexID(vertex));
  if (value == RN_UNKNOWN) return FALSE;
  for (int i = 0; i < mesh->VertexValence(vertex); i++) {
    R3MeshEdge *edge = mesh->EdgeOnVertex(vertex, i);
    R3MeshVertex *neighbor_vertex = mesh->VertexAcrossEdge(edge, vertex);
    RNScalar neighbor_value = property->VertexValue(mesh->VertexID(neighbor_vertex));
    if (neighbor_value == RN_UNKNOWN) continue;
    minimum &= (value + epsilon < neighbor_value) ? TRUE : FALSE;
    if ((extrema_type == PROPERTY_MINIMA) && !minimum) return FALSE;
    maximum &= (value - epsilon > neighbor_value) ? TRUE : FALSE;
    if ((extrema_type == PROPERTY_MAXIMA) && !maximum) return FALSE;
  }

  // Vertex value is extremum
  if (extrema_type == PROPERTY_MINIMA) return minimum;
  else if (extrema_type == PROPERTY_MAXIMA) return maximum;
  return minimum || maximum;
}



static RNArray<Point *> *
SelectPropertyExtrema(R3Mesh *mesh, R3MeshProperty *prp,
  int npoints, double min_spacing, int extrema_type)
{
  // Copy property 
  R3MeshProperty property(*prp);

  // Normalize property
  property.Normalize();

  // Select points from successively blurred property until number of extrema is npoints
  int max_blurs = 0;
  RNScalar epsilon = 0;
  RNScalar sigma = 0.01 * sqrt(mesh->Area());
  RNArray<R3MeshVertex *> selected_vertices;
  for (int blur = 0; blur <= max_blurs; blur++) { 
    // Find extremum vertices
    RNArray<R3MeshVertex *> extremum_vertices;
    for (int i = 0; i < mesh->NVertices(); i++) {
      R3MeshVertex *vertex = mesh->Vertex(i);
      if (!IsLocalExtremum(mesh, vertex, &property, extrema_type, epsilon)) continue; 
      extremum_vertices.Insert(vertex);
    }

    // Sort extremum vertices by property value
    for (int i = 0; i < extremum_vertices.NEntries(); i++) {
      for (int j = i+1; j < extremum_vertices.NEntries(); j++) {
        RNScalar value_i = property.VertexValue(extremum_vertices[i]);
        RNScalar value_j = property.VertexValue(extremum_vertices[j]);
        if (extrema_type == PROPERTY_MINIMA) { if (value_i > value_j) extremum_vertices.Swap(i, j); }
        else if (extrema_type == PROPERTY_MAXIMA) { if (value_i < value_j) extremum_vertices.Swap(i, j); }
        else { if (fabs(value_i) < fabs(value_j)) extremum_vertices.Swap(i, j); }
      }
    }

    // Select amongst extremum vertices taking into account min_spacing
    R3mesh_mark++;
    selected_vertices.Empty();
    for (int i = 0; i < extremum_vertices.NEntries(); i++) {
      R3MeshVertex *vertex = extremum_vertices.Kth(i);

      // Check if vertex is marked (i.e., too close to vertex previously selected)
      if (mesh->VertexMark(vertex) == R3mesh_mark) continue;

      // Check if vertex is disconnected completely
      if (mesh->VertexValence(vertex) == 0) continue;
    
      // Select vertex
      selected_vertices.Insert(vertex);

      // Mark nearby vertices
      mesh->SetVertexMark(vertex, R3mesh_mark);
      if (min_spacing > 0) {
        RNLength *distances = mesh->DijkstraDistances(vertex, min_spacing);
        for (int j = 0; j < mesh->NVertices(); j++) {
          R3MeshVertex *nearby_vertex = mesh->Vertex(j);
          if (mesh->VertexMark(nearby_vertex) != R3mesh_mark) {
            if (distances[j] <= min_spacing) {
              mesh->SetVertexMark(nearby_vertex, R3mesh_mark);
            }
          }
        } 
        delete [] distances;
      }
    }

    // Print debug statement
    if (print_debug) printf("  %d : %d %g %g\n", blur, selected_vertices.NEntries(), sigma, property.StandardDeviation());

    // Check if found correct number of points
    if ((npoints <= 0) || (selected_vertices.NEntries() <= npoints)) break;

    // Blur the property and iterate
    if (blur < max_blurs) {
      property.Blur(sigma);
      sigma *= 1.1;
    }
  }

  // For debugging purposes
  if (print_debug) property.Write("foo.val");

  // Allocate array of points
  RNArray<Point *> *points = new RNArray<Point *>();
  if (!points) {
    RNFail("Unable to allocate array of points\n");
    return NULL;
  }

  // Create points
  for (int i = 0; (i < selected_vertices.NEntries()) && (i < npoints); i++) {
    R3MeshVertex *vertex = selected_vertices[i];
    Point *point = new Point(mesh, vertex);
    points->Insert(point);
  }

  // Return points
  return points;
}



////////////////////////////////////////////////////////////////////////
// Scale-space sampling
////////////////////////////////////////////////////////////////////////

struct ScaleSpaceScore {
  int vertex_index;
  double score;
};
  

static int
CompareScaleSpaceScores(const void *s1, const void *s2)
{
  const ScaleSpaceScore *score1 = (const ScaleSpaceScore *) s1;
  const ScaleSpaceScore *score2 = (const ScaleSpaceScore *) s2;
  RNScalar delta = score1->score - score2->score;
  if (delta > 0) return -1;
  else if (delta < 0) return 1;
  else return 0;
}



static RNScalar 
VoteValue(R3Mesh *mesh, R3MeshVertex *vertex, R3MeshProperty *property, int extrema_type)
{
  // Get vertex value
  RNScalar value = property->VertexValue(vertex);
  if (value == RN_UNKNOWN) return 0;

  // Compute statistics of neighbors
  int count = 0;
  RNScalar sum_value = 0;
  RNScalar minimum_value = FLT_MAX;
  RNScalar maximum_value = -FLT_MAX;
  for (int i = 0; i < mesh->VertexValence(vertex); i++) {
    R3MeshEdge *edge = mesh->EdgeOnVertex(vertex, i);
    R3MeshVertex *neighbor_vertex = mesh->VertexAcrossEdge(edge, vertex);
    RNScalar neighbor_value = property->VertexValue(mesh->VertexID(neighbor_vertex));
    if (neighbor_value == RN_UNKNOWN) continue;
    if (neighbor_value < minimum_value) minimum_value = neighbor_value;
    if (neighbor_value > maximum_value) maximum_value = neighbor_value;
    sum_value += neighbor_value;
    count++;
  }

  // Just checking
  if (count == 0) return 0;

  // Return vote value
  RNScalar mean_value = sum_value / count;
  RNScalar laplacian = value - mean_value;
  if ((extrema_type == PROPERTY_MINIMA) && (value >= minimum_value)) return 0;
  if ((extrema_type == PROPERTY_MAXIMA) && (value <= maximum_value)) return 0;
  return fabs(laplacian);
}



static RNArray<Point *> *
SelectScaleSpaceExtrema(R3Mesh *mesh, R3MeshProperty *prp,
  int npoints, double min_spacing, int extrema_type)
{
  // Copy property 
  R3MeshProperty property(*prp);

  // Normalize property
  property.Normalize();

  // Allocate scores
  ScaleSpaceScore *scores = new ScaleSpaceScore [ mesh->NVertices() ];
  for (int i = 0; i < mesh->NVertices(); i++) {
    scores[i].vertex_index = i;
    scores[i].score = 0;
  }

  // Compute scale space scores
  int nblurs = 4;
  RNScalar weight = 1;
  RNScalar sigma = 0.01 * sqrt(mesh->Area());
  for (int blur = 0; blur < nblurs; blur++) { 
    // Compute strength of every vertex
    R3MeshProperty vote(property);
    for (int i = 0; i < mesh->NVertices(); i++) {
      R3MeshVertex *vertex = mesh->Vertex(i);
      RNScalar value = VoteValue(mesh, vertex, &property, extrema_type);
      vote.SetVertexValue(i, value); 
    }

    // Vote for vertex with its normalized strength value (augmented if local extremum)
    for (int i = 0; i < mesh->NVertices(); i++) {
      scores[i].score += weight * vote.VertexValue(i);
    }

    // Blur the property and iterate
    property.Blur(sigma);
    // weight *= 2;
    // sigma *= 2;
  }


  if (print_debug) {
    property.Write("property.val");
    R3MeshProperty foo(mesh, "ScaleSpaceScore");
    for (int i = 0; i < mesh->NVertices(); i++) foo.SetVertexValue(i, scores[i].score);
    foo.Write("scores.val");
  }

  // Sort vertices by score
  qsort(scores, mesh->NVertices(), sizeof(ScaleSpaceScore), CompareScaleSpaceScores);

  // Select vertices taking into account min_spacing
  R3mesh_mark++;
  RNArray<R3MeshVertex *> selected_vertices;
  for (int i = 0; i < mesh->NVertices(); i++) {
    R3MeshVertex *vertex = mesh->Vertex(scores[i].vertex_index);
    
    // Check if vertex is marked (i.e., too close to vertex previously selected)
    if (mesh->VertexMark(vertex) == R3mesh_mark) continue;

    // Check if vertex is disconnected completely
    if (mesh->VertexValence(vertex) == 0) continue;
    
    // Select vertex
    selected_vertices.Insert(vertex);
    
    // Mark nearby vertices
    mesh->SetVertexMark(vertex, R3mesh_mark);
    if (min_spacing > 0) {
      RNLength *distances = mesh->DijkstraDistances(vertex, min_spacing);
      for (int j = 0; j < mesh->NVertices(); j++) {
        R3MeshVertex *nearby_vertex = mesh->Vertex(j);
        if (mesh->VertexMark(nearby_vertex) != R3mesh_mark) {
          if (distances[j] <= min_spacing) {
            mesh->SetVertexMark(nearby_vertex, R3mesh_mark);
          }
        }
      } 
      delete [] distances;
    }

    // Check if found all points
    if ((npoints > 0) && (selected_vertices.NEntries() >= npoints)) break;
  }

  // Delete scores
  delete [] scores;

  // Allocate array of points
  RNArray<Point *> *points = new RNArray<Point *>();
  if (!points) {
    RNFail("Unable to allocate array of points\n");
    return NULL;
  }

  // Create points
  for (int i = 0; i < selected_vertices.NEntries(); i++) {
    R3MeshVertex *vertex = selected_vertices[i];
    Point *point = new Point(mesh, vertex);
    points->Insert(point);
  }

  // Return points
  return points;
}



////////////////////////////////////////////////////////////////////////
// Exterior point sampling
////////////////////////////////////////////////////////////////////////

static RNArray<Point *> *
SelectExteriorSurfacePoints(R3Mesh *mesh, int npoints, double min_spacing)
{
  // Allocate array of points
  RNArray<Point *> *points = new RNArray<Point *>();
  if (!points) {
    RNFail("Unable to allocate array of points\n");
    return NULL;
  }

  // Initialize grid
  RNArea area_per_point = mesh->Area() / npoints;
  RNScalar grid_spacing = sqrt(area_per_point);
  if (grid_spacing < min_spacing) grid_spacing = min_spacing;
  R3Grid exterior_grid(mesh->BBox(), grid_spacing, 5, 1024, 2);
  const RNScalar *exterior_grid_values = exterior_grid.GridValues();
  grid_spacing = exterior_grid.GridToWorldScaleFactor();

  // Initialize mesh search tree
  R3MeshSearchTree search_tree(mesh);

  // Initialze point search tree
  Point tmp; int position_offset = (unsigned char *) &(tmp.position) - (unsigned char *) &tmp;
  R3Kdtree<Point *> kdtree(exterior_grid.WorldBox(), position_offset);

  // Seed search with border voxels (assumed outside)
  int ix, iy, iz, grid_index;
  RNArray<const RNScalar *> stack;
  for (int i = 0; i < exterior_grid.XResolution(); i++) {
    for (int j = 0; j < exterior_grid.YResolution(); j++) {
      for (int k = 0; k < exterior_grid.ZResolution(); k++) {
        if ((i == 0) || (i == exterior_grid.XResolution()-1) ||
            (j == 0) || (j == exterior_grid.YResolution()-1) ||
            (k == 0) || (k == exterior_grid.ZResolution()-1)) {
          exterior_grid.IndicesToIndex(i, j, k, grid_index);
          stack.Insert(&exterior_grid_values[grid_index]);
          exterior_grid.SetGridValue(grid_index, 1.0);
        }
      }
    }
  }

  // Flood fill
  while (!stack.IsEmpty()) {
    const RNScalar *exterior_grid_valuep = stack.Tail(); stack.RemoveTail();
    int grid_index = exterior_grid_valuep - exterior_grid_values;
    assert((grid_index >= 0) && (grid_index < exterior_grid.NEntries()));
    assert(exterior_grid.GridValue(grid_index) > 0.5);
    exterior_grid.IndexToIndices(grid_index, ix, iy, iz);
    R3Point p0 = exterior_grid.WorldPosition(ix, iy, iz);
    for (int dz = -1; dz <= 1; dz++) {
      int gz = iz + dz;
      if ((gz < 0) || (gz >= exterior_grid.ZResolution())) continue;
      for (int dy = -1; dy <= 1; dy++) {
        int gy = iy + dy;
        if ((gy < 0) || (gy >= exterior_grid.YResolution())) continue;
        for (int dx = -1; dx <= 1; dx++) {
          int gx = ix + dx;
          if ((gx < 0) || (gx >= exterior_grid.XResolution())) continue;
          if (exterior_grid.GridValue(gx, gy, gz) > 0.5) continue;

#if 1       
          // Create point sample at front-most surface boundary separating cell from its neighbors
          RNBoolean blocked = FALSE;
          R3Point p1 = exterior_grid.WorldPosition(gx, gy, gz);
          R3Ray ray(p0, p1);
          R3Point midpoint = 0.5*(p0 + p1);
          RNArray<R3MeshIntersection *> hits;
          search_tree.FindAll(midpoint, hits, 0, grid_spacing);
          RNScalar best_t = FLT_MAX;
          R3Point best_position = R3zero_point;
          R3Vector best_normal = R3zero_vector;
          for (int k = 0; k < hits.NEntries(); k++) {
            R3MeshIntersection *hit = hits.Kth(k);
            R3Plane plane = mesh->FacePlane(hit->face);
            RNScalar d0 = R3SignedDistance(plane, p0);
            RNScalar d1 = R3SignedDistance(plane, p1);
            if (RNIsNegativeOrZero(d0*d1)) {
              if (mesh->Intersection(ray, hit->face, hit)) {
                if (hit->t < best_t) {
                  best_position = hit->point;
                  best_normal = plane.Normal();
                  if (d0 < 0) best_normal.Flip();
                  best_t = hit->t;
                }
              }
            }
            delete hit;
          }
          if (best_t < FLT_MAX) {
            if ((min_spacing == 0) || !kdtree.FindAny(best_position, 0, min_spacing)) {
              Point *point = new Point(best_position, best_normal);
              kdtree.InsertPoint(point);
              points->Insert(point);
            }
            blocked = TRUE;
          }
#else
          // FOR SOME REASON, THIS DOES NOT WORK
          // Create point sample at front-most surface boundary separating cell from its neighbors
          RNBoolean blocked = FALSE;
          R3Point p1 = exterior_grid.WorldPosition(gx, gy, gz);
          R3Ray ray(p0, p1);
          R3MeshIntersection hit;
          search_tree.FindIntersection(ray, hit, 0, grid_spacing);
          if (hit.type != R3_MESH_NULL_TYPE) {
            if ((min_spacing == 0) || !kdtree.FindAny(hit.point, 0, min_spacing)) {
              R3Point p = hit.point;
              R3Vector n = mesh->FaceNormal(hit.face);
              if (n.Dot(ray.Vector()) > 0) n.Flip();
              Point *point = new Point(p, n);
              kdtree.InsertPoint(point);
              points->Insert(point);
            }
            blocked = TRUE;
          }
#endif
          
          // Continue DFS 
          if (!blocked) {
            exterior_grid.IndicesToIndex(gx, gy, gz, grid_index);
            stack.Insert(&exterior_grid_values[grid_index]);
            exterior_grid.SetGridValue(grid_index, 1.0);
          }
        }
      }
    }
  }

  // exterior_grid.WriteFile("exterior.grd");

  // Return points
  return points;
}



////////////////////////////////////////////////////////////////////////
// Visible point sampling
////////////////////////////////////////////////////////////////////////

static RNArray<Point *> *
SelectVisibleSurfacePoints(R3Mesh *mesh, int npoints, double min_spacing)
{
  // Allocate array of points
  RNArray<Point *> *points = new RNArray<Point *>();
  if (!points) {
    RNFail("Unable to allocate array of points\n");
    return NULL;
  }

  // Initialze point search tree
  Point tmp; int position_offset = (unsigned char *) &(tmp.position) - (unsigned char *) &tmp;
  R3Kdtree<Point *> kdtree(mesh->BBox(), position_offset);

  // Get convenient variables
  // Need high resolution to avoid z-buffer quantization issues
  int resolution = 1024; 
  int nviews = npoints / (resolution * resolution);
  if (nviews < 128) nviews = 128;
  int npixels = resolution * resolution * nviews;
  RNScalar pixel_probability = (RNScalar) npoints / (RNScalar) npixels;
  R3Point c = mesh->BBox().Centroid();
  RNScalar r = mesh->BBox().DiagonalRadius();
  R2Box bbox2d(-r, -r, r, r);
  R2Grid face_image(bbox2d, r / resolution);
  R2Grid depth_image(face_image);

  // Insert points visible to random views
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

    // Generate points
    for (int iy = 0; iy < depth_image.YResolution(); iy++) {
      for (int ix = 0; ix < depth_image.XResolution(); ix++) {
        if (RNRandomScalar() > pixel_probability) continue;
        RNScalar face_value = face_image.GridValue(ix, iy);
        if (face_value < 0)  continue;
        RNScalar depth = depth_image.GridValue(ix, iy);
        if (depth == RN_INFINITY) continue;
        R2Point xy_position = depth_image.WorldPosition(ix, iy);
        R3Point position(xy_position.X(), xy_position.Y(), depth);
        position.InverseTransform(transformation);
        if ((min_spacing > 0) && kdtree.FindAny(position, 0, min_spacing)) continue;
        int face_index = (int) (face_value + 0.5);
        if ((face_index < 0) || (face_index >= mesh->NFaces())) continue;
        R3MeshFace *face = mesh->Face(face_index);
        R3Vector normal = mesh->FaceNormal(face);
        R3Vector n = normal; n.Transform(transformation);
        if (n.Z() > 0) normal.Flip();
        Point *point = new Point(position, normal);
        if (min_spacing > 0) kdtree.InsertPoint(point);
        points->Insert(point);
      }
    }
  }

  // Return points
  return points;
}



////////////////////////////////////////////////////////////////////////
// Near surface sampling
////////////////////////////////////////////////////////////////////////

static RNArray<Point *> *
SelectNearSurfacePoints(R3Mesh *mesh, R3MeshProperty *property, int npoints)
{
  // Allocate array of points
  RNArray<Point *> *points = new RNArray<Point *>();
  if (!points) {
    RNFail("Unable to allocate array of points\n");
    return NULL;
  }

  // Create kdtree
  R3MeshSearchTree kdtree(mesh);

  // Compute face sampling weights
  RNScalar total_weight = ComputeFaceSamplingWeights(mesh, property);
  if (total_weight == 0) return NULL;

  // Sample points randomly from mesh surface
  for (int i = 0; i < mesh->NFaces(); i++) {
    R3MeshFace *face = mesh->Face(i);
    RNScalar face_weight = mesh->FaceValue(face);

    // Determine number of points on face (twice needed)
    RNScalar ideal_face_npoints = 1.1 * 0.5 * npoints * face_weight / total_weight;
    int face_npoints = (int) ideal_face_npoints;
    RNScalar remainder = ideal_face_npoints - face_npoints;
    if (remainder > RNRandomScalar()) face_npoints++;

    // Sample points at random position on face
    for (int j = 0; j < face_npoints; j++) {
      R3Point position = mesh->RandomPointOnFace(face);

      // Get ray from point in normal direction
      R3Vector normal = R3zero_vector;
      R3Point barycentrics = mesh->FaceBarycentric(face, position);
      R3MeshVertex *v0 = mesh->VertexOnFace(face, 0);
      R3MeshVertex *v1 = mesh->VertexOnFace(face, 1);
      R3MeshVertex *v2 = mesh->VertexOnFace(face, 2);
      normal += barycentrics[0] * mesh->VertexNormal(v0);
      normal += barycentrics[1] * mesh->VertexNormal(v1);
      normal += barycentrics[2] * mesh->VertexNormal(v2);
      normal.Normalize();
  
      // Find ray intersection on outside
      R3MeshIntersection intersectionA;
      R3Ray rayA(position, normal);
      kdtree.FindIntersection(rayA, intersectionA, RN_EPSILON, max_distance);
      if (intersectionA.type == R3_MESH_NULL_TYPE) {
        intersectionA.t = max_distance;
      }

      // Find ray intersection on inside
      R3MeshIntersection intersectionB;
      R3Ray rayB(position, -normal);
      kdtree.FindIntersection(rayB, intersectionB, RN_EPSILON, max_distance);
      if (intersectionB.type == R3_MESH_NULL_TYPE) {
        intersectionB.t = max_distance;
      }

      // Find lesser of two ray intersection distances
      RNScalar t = intersectionA.t;
      if (intersectionB.t < t) t = intersectionB.t;
      if (t == FLT_MAX) continue;  // should not happen for watertight mesh or if max_distance is set

      // Adjust t based on max distance
      if (t > max_distance) t = max_distance;

      // Adjust t to bias towards smaller values
      t *= pow(RNRandomScalar(), near_surface_bias_exponent);

      // Create point on outside
      R3Point positionA = position + t * normal;
      if (bbox.IsEmpty() || R3Contains(bbox, positionA)) {
        // Compute sdf
        RNScalar sdfA = 1;
        if (!binary_sdf) {
          R3MeshIntersection closest;
          kdtree.FindClosest(positionA, closest, 0, max_distance);
          sdfA = (closest.type != R3_MESH_NULL_TYPE) ? closest.t : t;
        }

        // Create point
        Point *pointA = new Point();
        pointA->face = face;
        pointA->position = positionA;
        pointA->normal = normal;
        pointA->sdf = sdfA;
        points->Insert(pointA);
      }

      // Create point on inside
      R3Point positionB = position - t * normal;
      if (bbox.IsEmpty() || R3Contains(bbox, positionB)) {
        // Compute sdf
        RNScalar sdfB = -1;
        if (!binary_sdf) {
          R3MeshIntersection closest;
          kdtree.FindClosest(positionB, closest, 0, max_distance);
          sdfB = (closest.type != R3_MESH_NULL_TYPE) ? -closest.t : -t;
        }

        R3MeshIntersection closest;
        kdtree.FindClosest(positionB, closest, 0, max_distance);
        Point *pointB = new Point();
        pointB->face = face;
        pointB->position = positionB;
        pointB->normal = normal;
        pointB->sdf = sdfB;
        points->Insert(pointB);
      }
    }
  }

  // Return points
  return points;
}



////////////////////////////////////////////////////////////////////////
// Uniform bbox point sampling
////////////////////////////////////////////////////////////////////////

static RNArray<Point *> *
SelectPointsUniformInBBox(R3Mesh *mesh, int npoints)
{
  // Get bounding box
  R3Box b = bbox;
  if (b.IsEmpty()) {
    b = mesh->BBox();
    b.Inflate(1.1);
  }
  
  // Create signed distance grid
  R3Grid *signed_distance_grid = CreateSignedDistanceGrid(mesh, b);
  if (!signed_distance_grid) {
    RNFail("Unable to create signed distance grid\n");
    return NULL;
  }

  // Allocate array of points
  RNArray<Point *> *points = new RNArray<Point *>();
  if (!points) {
    RNFail("Unable to allocate array of points\n");
    return NULL;
  }
  
  // Create kdtree
  R3MeshSearchTree *kdtree = NULL;
  if (!binary_sdf) {
    kdtree = new R3MeshSearchTree(mesh);
    assert(kdtree);
  }

  // Sample points randomly from mesh bbox interior
  for (int i = 0; i < npoints; i++) {
    // Sample location
    RNScalar x = RNRandomScalar() * b.XLength() + b.XMin();
    RNScalar y = RNRandomScalar() * b.YLength() + b.YMin();
    RNScalar z = RNRandomScalar() * b.ZLength() + b.ZMin();
    R3Point position(x, y, z);

    // Compute sign 
    RNScalar sdf = signed_distance_grid->WorldValue(position);
    RNScalar sign = (sdf < 0) ? -1 : 1;
    
    // Compute distance to closest surface
    if (kdtree) {
      R3MeshIntersection closest;
      kdtree->FindClosest(position, closest);
      sdf = sign * closest.t;
    }

    // Create point
    Point *point = new Point();
    point->position = position;
    point->normal = R3zero_vector;
    point->sdf = sdf;
    points->Insert(point);
  }

  // Delete kdtree
  if (kdtree) delete kdtree;

  // Delete signed distance grid
  delete signed_distance_grid;

  // Return points
  return points;
}



////////////////////////////////////////////////////////////////////////
// Select one point at center of mass
////////////////////////////////////////////////////////////////////////

static RNArray<Point *> *
SelectCenterOfMass(R3Mesh *mesh)
{
  // Allocate array of points
  RNArray<Point *> *points = new RNArray<Point *>();
  if (!points) {
    RNFail("Unable to allocate array of points\n");
    return NULL;
  }

  // Compute mesh center of mass
  R3Point mesh_centroid = mesh->Centroid();

  // Insert mesh centroid
  Point *point = new Point(mesh_centroid, R3zero_vector);
  points->Insert(point);

  // Return points
  return points;
}



////////////////////////////////////////////////////////////////////////
// Select one point at vertex closest to center of mass
////////////////////////////////////////////////////////////////////////

static RNArray<Point *> *
SelectVertexClosestToCenterOfMass(R3Mesh *mesh)
{
  // Allocate array of points
  RNArray<Point *> *points = new RNArray<Point *>();
  if (!points) {
    RNFail("Unable to allocate array of points\n");
    return NULL;
  }

  // Compute mesh center of mass
  R3Point mesh_centroid = mesh->Centroid();

  // Find closest vertex 
  RNScalar closest_d = FLT_MAX;
  R3MeshVertex *closest_vertex = NULL;
  for (int i = 0; i < mesh->NVertices(); i++) {
    R3MeshVertex *vertex = mesh->Vertex(i);
    const R3Point& position = mesh->VertexPosition(vertex);
    RNScalar d = R3Distance(position, mesh_centroid);
    if (d < closest_d) {
      closest_vertex = vertex;
      closest_d = d;
    }
  }

  // Insert closest vertex
  assert(closest_vertex);
  Point *point = new Point(mesh, closest_vertex);
  points->Insert(point);

  // Return points
  return points;
}



////////////////////////////////////////////////////////////////////////
// Top-level point sampling function
////////////////////////////////////////////////////////////////////////

static RNArray<Point *> *
SelectPoints(R3Mesh *mesh, R3MeshProperty *property, int selection_method)
{
  // Start statistics
  RNTime start_time;
  start_time.Read();

  // Update min spacing 
  if ((min_spacing < 0) && (min_normalized_spacing > 0)) {
    min_spacing = min_normalized_spacing * sqrt(mesh->Area());
  }

  // Update number of points
  if (num_points <= 0) {
    if (min_spacing > 0) num_points = mesh->NVertices();
    else num_points = 1024;
  }

  // Update min spacing
  if ((num_points > 0) && (min_spacing <= 0) && (min_relative_spacing > 0)) {
    min_spacing = min_relative_spacing * 2 * sqrt(mesh->Area() / (RN_PI * num_points));
  }

  // Select points using requested selection method 
  RNArray<Point *> *points = NULL;
  if (selection_method == RANDOM_VERTICES) 
    points = SelectVertexPoints(mesh, num_points, min_spacing);
  else if (selection_method == RANDOM_SURFACE_POINTS) 
    points= SelectRandomSurfacePoints(mesh, property, num_points);
  else if (selection_method == EXTERIOR_SURFACE_POINTS) 
    points= SelectExteriorSurfacePoints(mesh, num_points, min_spacing);
  else if (selection_method == NEAR_SURFACE_POINTS) 
    points= SelectNearSurfacePoints(mesh, property, num_points);
  else if (selection_method == VISIBLE_SURFACE_POINTS) 
    points= SelectVisibleSurfacePoints(mesh, num_points, min_spacing);
  else if (selection_method == ITERATIVE_FURTHEST_VERTEX) 
    points= SelectFurthestPoints(mesh, num_points, min_spacing);
  else if (selection_method == PROPERTY_MINIMA) 
    points= SelectPropertyExtrema(mesh, property, num_points, min_spacing, PROPERTY_MINIMA);
  else if (selection_method == PROPERTY_MAXIMA) 
    points= SelectPropertyExtrema(mesh, property, num_points, min_spacing, PROPERTY_MAXIMA);
  else if (selection_method == PROPERTY_EXTREMA)
    points= SelectPropertyExtrema(mesh, property, num_points, min_spacing, PROPERTY_EXTREMA);
  else if (selection_method == SCALE_SPACE_MINIMA) 
    points= SelectScaleSpaceExtrema(mesh, property, num_points, min_spacing, PROPERTY_MINIMA);
  else if (selection_method == SCALE_SPACE_MAXIMA) 
    points= SelectScaleSpaceExtrema(mesh, property, num_points, min_spacing, PROPERTY_MAXIMA);
  else if (selection_method == SCALE_SPACE_EXTREMA)
    points= SelectScaleSpaceExtrema(mesh, property, num_points, min_spacing, PROPERTY_EXTREMA);
  else if (selection_method == UNIFORM_IN_BBOX) 
    points= SelectPointsUniformInBBox(mesh, num_points);
  else if (selection_method == CENTER_OF_MASS) 
    points= SelectCenterOfMass(mesh);
  else if (selection_method == VERTEX_CLOSEST_TO_CENTER_OF_MASS) 
    points= SelectVertexClosestToCenterOfMass(mesh);

  // Check points
  if (!points) {
    RNFail("Error selecting points\n");
    return 0;
  }

  // Check if too few points
  if (min_points > 0) {
    while (points->NEntries() < min_points) {
      int face_index = (int) (RNRandomScalar() * mesh->NFaces());
      R3MeshFace *face = mesh->Face(face_index);
      R3Point position = mesh->RandomPointOnFace(face);
      const R3Vector& normal = mesh->FaceNormal(face);
      Point *point = new Point(position, normal);
      points->Insert(point);
    }
  }

  // Check if too many points
  if (max_points > 0) {
    while (points->NEntries() > max_points) {
      int k = RNRandomScalar() * points->NEntries();
      points->Swap(k, points->NEntries()-1);
      Point *point = points->Tail();
      points->RemoveTail();
      delete point;
    }
  }

  // Print message
  if (print_verbose) {
    fprintf(stdout, "Generated points ...\n");
    printf("  Time = %.2f seconds\n", start_time.Elapsed());
    printf("  # Points = %d\n", points->NEntries());
    printf("  Selection method = %d\n", selection_method);
    if (min_points > 0) printf("  Minimum # points = %d\n", min_points);
    if (max_points > 0) printf("  Maximum # points = %d\n", max_points);
    printf("  Requested # points = %d\n", num_points);
    printf("  Requested min spacing = %g\n", min_spacing);
    fflush(stdout);
  }

  // Return array of points
  return points;
}



////////////////////////////////////////////////////////////////////////
// Program argument parsing
////////////////////////////////////////////////////////////////////////

static int
ParseArgs (int argc, char **argv)
{
  // Parse arguments
  argc--; argv++;
  while (argc > 0) {
    if ((*argv)[0] == '-') {
      if (!strcmp(*argv, "-num_points")) { argv++; argc--; min_points = max_points = num_points = atoi(*argv); }
      else if (!strcmp(*argv, "-npoints")) { argv++; argc--; num_points = atoi(*argv); }
      else if (!strcmp(*argv, "-min_points")) { argv++; argc--; min_points = atoi(*argv); }
      else if (!strcmp(*argv, "-max_points")) { argv++; argc--; max_points = atoi(*argv); }
      else if (!strcmp(*argv, "-min_spacing")) { argc--; argv++; min_spacing = atof(*argv); }
      else if (!strcmp(*argv, "-min_normalized_spacing")) { argc--; argv++; min_normalized_spacing = atof(*argv); }
      else if (!strcmp(*argv, "-min_relative_spacing")) { argc--; argv++; min_relative_spacing = atof(*argv); }
      else if (!strcmp(*argv, "-max_distance")) { argc--; argv++; max_distance = atof(*argv); }
      else if (!strcmp(*argv, "-curvature_exponent")) { argc--; argv++; curvature_exponent = atof(*argv); }
      else if (!strcmp(*argv, "-curvature_max")) { argc--; argv++; curvature_max = atof(*argv); }
      else if (!strcmp(*argv, "-near_surface_bias_exponent")) { argc--; argv++; near_surface_bias_exponent = atof(*argv); }
      else if (!strcmp(*argv, "-property")) { argc--; argv++;  property_name = *argv; }
      else if (!strcmp(*argv, "-selection_method")) { argc--; argv++; selection_method = atoi(*argv); }
      else if (!strcmp(*argv, "-random_surface_points")) { selection_method = RANDOM_SURFACE_POINTS; }
      else if (!strcmp(*argv, "-exterior_surface_points")) { selection_method = EXTERIOR_SURFACE_POINTS; }
      else if (!strcmp(*argv, "-visible_surface_points")) { selection_method = VISIBLE_SURFACE_POINTS; }
      else if (!strcmp(*argv, "-near_surface")) { selection_method = NEAR_SURFACE_POINTS; }
      else if (!strcmp(*argv, "-near_surface_points")) { selection_method = NEAR_SURFACE_POINTS; }
      else if (!strcmp(*argv, "-random_vertices")) { selection_method = RANDOM_VERTICES; }
      else if (!strcmp(*argv, "-iterative_furthest_vertex")) { selection_method = ITERATIVE_FURTHEST_VERTEX; }
      else if (!strcmp(*argv, "-center_of_mass")) { selection_method = CENTER_OF_MASS; }
      else if (!strcmp(*argv, "-vertex_closest_to_center_of_mass")) { selection_method = VERTEX_CLOSEST_TO_CENTER_OF_MASS; }
      else if (!strcmp(*argv, "-property_minima")) { selection_method = PROPERTY_MINIMA; }
      else if (!strcmp(*argv, "-property_maxima")) { selection_method = PROPERTY_MAXIMA; }
      else if (!strcmp(*argv, "-property_extrema")) { selection_method = PROPERTY_EXTREMA; }
      else if (!strcmp(*argv, "-scale_space_minima")) { selection_method = SCALE_SPACE_MINIMA; }
      else if (!strcmp(*argv, "-scale_space_maxima")) { selection_method = SCALE_SPACE_MAXIMA; }
      else if (!strcmp(*argv, "-scale_space_extrema")) { selection_method = SCALE_SPACE_EXTREMA; }
      else if (!strcmp(*argv, "-uniform_in_bbox")) { selection_method = UNIFORM_IN_BBOX; }
      else if (!strcmp(*argv, "-binary_sdf")) { binary_sdf = TRUE; }
      else if (!strcmp(*argv, "-v")) { print_verbose = TRUE; }
      else if (!strcmp(*argv, "-debug")) { print_debug = TRUE; }
      else if (!strcmp(*argv, "-bbox")) {
        argc--; argv++; bbox[0][0] = atof(*argv);
        argc--; argv++; bbox[0][1] = atof(*argv);
        argc--; argv++; bbox[0][2] = atof(*argv);
        argc--; argv++; bbox[1][0] = atof(*argv);
        argc--; argv++; bbox[1][1] = atof(*argv);
        argc--; argv++; bbox[1][2] = atof(*argv);
      }
      else { RNFail("Invalid program argument: %s\n", *argv); return 0; }
      argv++; argc--;
    }
    else {
      if (!mesh_name) mesh_name = *argv;
      else if (!points_name) points_name = *argv;
      else RNAbort("Invalid program argument: %s", *argv);
      argv++; argc--;
    }
  }

  // Check mesh and points names
  if (!mesh_name || !points_name) {
    // Print usage statement
    RNFail("Usage: msh2pts <options> input_mesh output_points\n\n");
    return FALSE;
  }

  // Resolve binary sdf flag
  if (!strstr(points_name, ".sdf")) binary_sdf = TRUE;

  // Return success
  return TRUE;
}



////////////////////////////////////////////////////////////////////////

int
main(int argc, char **argv)
{
  // Start statistics
  RNTime start_time;
  start_time.Read();

  // Parse args
  if(!ParseArgs(argc, argv)) exit(-1);

  // Read mesh
  R3Mesh *mesh = ReadMesh(mesh_name);
  if (!mesh) exit(-1);

  // Read property 
  R3MeshProperty *property = NULL;
  if (property_name) {
    property = ReadProperty(mesh, property_name);
    if (!property) exit(-1);
  }

  // Select points
  RNArray<Point *> *points = SelectPoints(mesh, property, selection_method);
  if (!points) exit(-1);

  // Write points
  if (!WritePoints(mesh, *points, points_name)) exit(-1);

  // Print message
  if (print_verbose) {
    fprintf(stdout, "Finished in %6.3f seconds.\n", start_time.Elapsed());
    fflush(stdout);
  }

  // Return success
  return 0;
}








