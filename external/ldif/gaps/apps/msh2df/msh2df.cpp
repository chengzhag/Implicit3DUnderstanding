// Source file for the mesh to distance function program



// Include files 

namespace gaps {}
using namespace gaps;
#include "R3Shapes/R3Shapes.h"



// Program variables

static char *input_mesh_filename = NULL;
static char *output_grid_filename = NULL;
static char *output_mesh_filename = NULL;
static char *output_points_filename = NULL;
static double target_grid_spacing = 0.01;  // in world units
static R3Box grid_bbox(FLT_MAX, FLT_MAX, FLT_MAX, -FLT_MAX, -FLT_MAX, -FLT_MAX);
static double grid_border = 3; // in grid units
static int grid_max_resolution = 512; // in grid units
static int refinement_radius = 9; // in grid units
static double truncation_distance = RN_INFINITY; // in world units
static int estimate_sign = FALSE;
static int sign_estimation_method = 1; // flood fill
static RNBoolean input_is_manifold = 0;
static RNBoolean input_is_range_scan = 0;
static R3Point scan_viewpoint(0,0,0);
static int print_verbose = 0;
static int print_debug = 0;



////////////////////////////////////////////////////////////////////////
// Type definitions
////////////////////////////////////////////////////////////////////////

struct Point {
  Point(void) : position(0,0,0), normal(0,0,0), radius(0) {};
  Point(const R3Point& position, const R3Vector& normal, RNLength radius)
    : position(position), normal(normal), radius(radius) {};
  R3Point position;
  R3Vector normal;
  RNLength radius;
};
  

  
////////////////////////////////////////////////////////////////////////
// Input/output
////////////////////////////////////////////////////////////////////////

static R3Mesh *
ReadMesh(char *filename)
{
  // Start statistics
  RNTime start_time;
  start_time.Read();

  // Allocate mesh
  R3Mesh *mesh = new R3Mesh();
  if (!mesh) {
    RNFail("Unable to allocate mesh for %s\n", filename);
    return NULL;
  }

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



static int 
WriteGrid(R3Grid *grid, const char *filename)
{
  // Start statistics
  RNTime start_time;
  start_time.Read();

  // Write grid
  if (!grid->WriteFile(filename)) return 0;

  // Print statistics
  if (print_verbose) {
    printf("Wrote grid ...\n");
    printf("  Time = %.2f seconds\n", start_time.Elapsed());
    printf("  # MBytes = %g\n", (double) (grid->NEntries() * sizeof(RNScalar)) / (1024.0 * 1024.0));
    fflush(stdout);
  }

  // Return success
  return 1;
}



static int
WriteMesh(R3Mesh *mesh, const char *filename)
{
  // Start statistics
  RNTime start_time;
  start_time.Read();

  // Write mesh
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
WritePoints(const RNArray<Point *>& points, const char *filename)
{
  // Start statistics
  RNTime start_time;
  start_time.Read();

  // Open file
  FILE *fp = fopen(output_points_filename, "wb");
  if (!fp) {
    RNFail("Unable to open points file %s\n", filename);
    return 0;
  }
  
  // Write points 
  float coordinates[6];
  for (int i = 0; i < points.NEntries(); i++) {
    Point *point = points.Kth(i);
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
  fclose(fp);
  
  // Print statistics
  if (print_verbose) {
    printf("Wrote points to %s ...\n", filename);
    printf("  Time = %.2f seconds\n", start_time.Elapsed());
    printf("  # Points = %d\n", points.NEntries());
    fflush(stdout);
  }

  // Return success
  return 1;
}



////////////////////////////////////////////////////////////////////////
// Distance estimation
////////////////////////////////////////////////////////////////////////

static int
EstimateDistanceUsingGrid(R3Grid *grid, R3Mesh *mesh)
{
  // Start statistics
  RNTime start_time;
  start_time.Read();

  // Rasterize each triangle into grid
  grid->Clear(0);
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
  
  // Print statistics
  if (print_verbose) {
    printf("Estimated distance using grid ...\n");
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

  // Return success
  return 1;
}



static int
RefineDistanceUsingKdtree(R3Grid *grid, R3Mesh *mesh, int refinement_radius)
{
  // Start statistics
  RNTime start_time;
  start_time.Read();

  // Create mesh search tree
  R3MeshSearchTree search_tree(mesh);

  // Compute the precise distance at every grid cell
  RNLength grid_spacing = grid->GridToWorldScaleFactor();
  RNScalar max_distance = refinement_radius * grid_spacing;
  if (max_distance > mesh->BBox().DiagonalLength()) max_distance = mesh->BBox().DiagonalLength();
  if (max_distance > truncation_distance) max_distance = truncation_distance;
  for (int iz = 0; iz < grid->ZResolution(); iz++) {
    for (int iy = 0; iy < grid->YResolution(); iy++) {
      for (int ix = 0; ix < grid->XResolution(); ix++) {
        R3Point world_position = grid->WorldPosition(ix, iy, iz);
        RNScalar grid_distance = grid->GridValue(ix, iy, iz);
        RNScalar sign = (grid_distance >= 0) ? 1 : -1;
        grid_distance = fabs(grid_distance);
        if (grid_distance > max_distance) continue;
        RNScalar distance = (grid_distance > grid_spacing) ? grid_distance : grid_spacing;
        while (distance < max_distance) {
          distance *= 1.25;
          R3MeshIntersection closest;
          search_tree.FindClosest(world_position, closest, 0, distance);
          if (closest.type != R3_MESH_NULL_TYPE) {
            grid->SetGridValue(ix, iy, iz, sign * closest.t);
            break;
          }
        }
      }
    }
  }
  
  // Print statistics
  if (print_verbose) {
    printf("Refined distance using kdtree ...\n");
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
    printf("  Refinement grid radius = %d\n", refinement_radius);
    printf("  Refinement world radius = %g\n", max_distance);
    fflush(stdout);
  }

  // Return success
  return 1;
}



////////////////////////////////////////////////////////////////////////
// Sign estimation
// -1=interior, 1=exterior
////////////////////////////////////////////////////////////////////////

static int
EstimateSignUsingNormals(R3Grid *grid, R3Mesh *mesh)
{
  // Initialize the sign grid
  R3Grid sign_grid(*grid);
  sign_grid.Clear(0);

  // Get convenient variables
  RNLength grid_spacing = grid->GridToWorldScaleFactor();

  // Rasterize negative sides of faces
  for (int i = 0; i < mesh->NFaces(); i++) {
    R3MeshFace *face = mesh->Face(i);
    R3Vector n = grid_spacing * mesh->FaceNormal(face);
    R3MeshVertex *v0 = mesh->VertexOnFace(face, 0);
    R3MeshVertex *v1 = mesh->VertexOnFace(face, 1);
    R3MeshVertex *v2 = mesh->VertexOnFace(face, 2);
    R3Point p0 = mesh->VertexPosition(v0);
    R3Point p1 = mesh->VertexPosition(v1);
    R3Point p2 = mesh->VertexPosition(v2);
    sign_grid.RasterizeWorldTriangle(p0 - n, p1 - n, p2 - n, -1.0, R3_GRID_REPLACE_OPERATION);
  }

  // Rasterize positive sides of faces (overwrite negatives if collisions)
  for (int i = 0; i < mesh->NFaces(); i++) {
    R3MeshFace *face = mesh->Face(i);
    R3Vector n = grid_spacing * mesh->FaceNormal(face);
    R3MeshVertex *v0 = mesh->VertexOnFace(face, 0);
    R3MeshVertex *v1 = mesh->VertexOnFace(face, 1);
    R3MeshVertex *v2 = mesh->VertexOnFace(face, 2);
    R3Point p0 = mesh->VertexPosition(v0);
    R3Point p1 = mesh->VertexPosition(v1);
    R3Point p2 = mesh->VertexPosition(v2);
    R3Vector n0 = grid_spacing * mesh->VertexNormal(v0);
    R3Vector n1 = grid_spacing * mesh->VertexNormal(v1);
    R3Vector n2 = grid_spacing * mesh->VertexNormal(v2);
    sign_grid.RasterizeWorldTriangle(p0, p1, p2, 1.0, R3_GRID_REPLACE_OPERATION);
    sign_grid.RasterizeWorldTriangle(p0 + n, p1 + n, p2 + n, 1.0, R3_GRID_REPLACE_OPERATION);
    sign_grid.RasterizeWorldTriangle(p0 + 0.5*n, p1 + 0.5*n, p2 + 0.5*n, 1.0, R3_GRID_REPLACE_OPERATION);
    sign_grid.RasterizeWorldTriangle(p0 + 1.5*n, p1 + 1.5*n, p2 + 1.5*n, 1.0, R3_GRID_REPLACE_OPERATION);
    sign_grid.RasterizeWorldTriangle(p0 + n0, p1 + n1, p2 + n2, 1.0, R3_GRID_REPLACE_OPERATION);
    sign_grid.RasterizeWorldTriangle(p0 + 0.5*n0, p1 + 0.5*n1, p2 + 0.5*n2, 1.0, R3_GRID_REPLACE_OPERATION);
    sign_grid.RasterizeWorldTriangle(p0 + 1.5*n0, p1 + 1.5*n1, p2 + 1.5*n2, 1.0, R3_GRID_REPLACE_OPERATION);
  }

  // Flood fill signs
  sign_grid.Voronoi();

  // Apply signs
  sign_grid.Threshold(0, -1, 1);
  grid->Multiply(sign_grid);

#if 0
  // Write some debug info
  if (print_debug) sign_grid.WriteFile("sign.grd");
#endif
  
  // Return success
  return 1;
}



static int
EstimateSignUsingFloodFill(R3Grid *grid, R3Mesh *mesh)
{
  // Initialize the sign grid
  R3Grid sign_grid(*grid);
  sign_grid.Clear(0);

  // Get convenient variables
  RNArray<const RNScalar *> stack;
  const RNScalar *sign_grid_values = sign_grid.GridValues();
  RNLength grid_spacing = grid->GridToWorldScaleFactor();
  RNScalar max_grid_step = sqrt(3)*grid_spacing;
  int ix, iy, iz, grid_index;

  // Initialize mesh search tree
  R3MeshSearchTree search_tree(mesh);
  
  // Seed search with border voxels (assumed outside)
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
    R3Point p0 = sign_grid.WorldPosition(ix, iy, iz);
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

          // Check if could step through surface
          RNScalar d = grid->GridValue(ix, iy, iz);
          if (d < max_grid_step) {
            RNScalar grid_step = sqrt((dx*dx + dy*dy + dz*dz) * grid_spacing*grid_spacing);
            if (d < grid_step) {
              // Check if stepping through a surface
              RNBoolean blocked = FALSE;
              R3Point p1 = sign_grid.WorldPosition(gx, gy, gz);
              R3Point midpoint = 0.5*(p0 + p1);
              RNArray<R3MeshIntersection *> hits;
              search_tree.FindAll(midpoint, hits, 0, 0.5*max_grid_step);
              for (int k = 0; k < hits.NEntries(); k++) {
                R3MeshIntersection *hit = hits.Kth(k);
                R3Plane plane = mesh->FacePlane(hit->face);
                delete hit;
                RNScalar d0 = R3SignedDistance(plane, p0);
                RNScalar d1 = R3SignedDistance(plane, p1);
                if (RNIsNegativeOrZero(d0*d1)) {
                  for (int k1 = k+1; k1 < hits.NEntries(); k1++) delete hits[k1];
                  blocked = TRUE;
                  break;
                }
              }
              if (blocked) continue;
            }
          }
             
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

#if 1
  // Write some debug info
  if (print_debug) sign_grid.WriteFile("sign.grd");
#endif
  
  // Return success
  return 1;
}



static int
EstimateSign(R3Grid *grid, R3Mesh *mesh)
{
  // Start statistics
  RNTime start_time;
  start_time.Read();


  // Check sign estimation method
  if (sign_estimation_method == 0) {
    if (!EstimateSignUsingNormals(grid, mesh)) return 0;
  }
  else if (sign_estimation_method == 1) {
    if (!EstimateSignUsingFloodFill(grid, mesh)) return 0;
  }
  else {
    RNFail("Unrecognized sign estimation method: %d\n", sign_estimation_method);
    return 0;
  }

  // Print statistics
  if (print_verbose) {
    printf("Estimated sign ...\n");
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

  // Return success
  return 1;
}



////////////////////////////////////////////////////////////////////////
// Point sampling
////////////////////////////////////////////////////////////////////////

static int
CreatePoints(R3Grid *grid, R3Mesh *mesh, RNArray<Point *>& points)
{
  // Start statistics
  RNTime start_time;
  start_time.Read();

  // Create mesh search tree
  R3MeshSearchTree search_tree(mesh);

  // Create point samples on boundary
  RNLength grid_spacing = grid->GridToWorldScaleFactor();
  RNLength max_distance = refinement_radius * grid_spacing;
  if (max_distance < grid_spacing) max_distance = grid_spacing;
  for (int iz = 0; iz < grid->ZResolution(); iz++) {
    for (int iy = 0; iy < grid->YResolution(); iy++) {
      for (int ix = 0; ix < grid->XResolution(); ix++) {

        // Check if in swath on positive side
        RNScalar grid_distance = grid->GridValue(ix, iy, iz);
        if (grid_distance < 0) continue;
        if (grid_distance > max_distance) continue;

        // Create point samples on boundary
        RNArray<R3MeshIntersection *> hits;
        R3Point world_position = grid->WorldPosition(ix, iy, iz);
        search_tree.FindAll(world_position, hits, 0, grid_distance + RN_EPSILON);
        for (int k = 0; k < hits.NEntries(); k++) {
          R3MeshIntersection *hit = hits.Kth(k);
          Point *point = new Point(hit->point, mesh->FaceNormal(hit->face), grid_spacing);
          points.Insert(point);
          delete hit;
        }
      }
    }
  }
    
  // Create kdtree for computing radii
  Point tmp; int position_offset = (unsigned char *) &(tmp.position) - (unsigned char *) &tmp;
  R3Kdtree<Point *> kdtree(points, position_offset);
  
  // Compute radius for each point
  int count = 0;
  RNScalar sum = 0;
  const int max_neighbors = 6;
  RNScalar neighbor_distances[max_neighbors];
  for (int i = 0; i < points.NEntries(); i++) {
    Point *point = points.Kth(i);
    RNArray<Point *> neighbors;
    point->radius = 0.5 * grid_spacing;
    if (kdtree.FindClosest(point->position, 0, 3*grid_spacing, max_neighbors, neighbors, neighbor_distances)) {
      if (!neighbors.IsEmpty()) {
        RNScalar radius = 0.5 * neighbor_distances[neighbors.NEntries()-1];
        if (radius < point->radius) {
          point->radius = radius;
          sum += radius;
          count++;
        }
      }
    }
  }

  // Print statistics
  if (print_verbose) {
    printf("Sampled points on boundary ...\n");
    printf("  Time = %.2f seconds\n", start_time.Elapsed());
    printf("  # Points = %d\n", points.NEntries());
    fflush(stdout);
  }

  // Return success
  return 1;
}



static int
DeletePoints(const RNArray<Point *>& points)
{
  // Delete point samples
  for (int i = 0;i < points.NEntries(); i++) delete points[i];

  // Return success
  return 1;
}



////////////////////////////////////////////////////////////////////////
// Signed distance update (recompute with respect to enclosing surface)
////////////////////////////////////////////////////////////////////////

static int
UpdateSignedDistanceUsingGrid(R3Grid *grid, R3Mesh *mesh)
{
  // Start statistics
  RNTime start_time;
  start_time.Read();

  // Recompute negative side
  R3Grid distance_grid(*grid);
  distance_grid.Threshold(0, 0, R3_GRID_KEEP_VALUE);
  distance_grid.SquaredDistanceTransform();
  distance_grid.Sqrt();
  distance_grid.Multiply(grid->GridToWorldScaleFactor());
  distance_grid.Negate();

  // Keep positive side, replace negative side
  for (int i = 0; i < grid->NEntries(); i++) {
    if (grid->GridValue(i) >= 0) continue;
    grid->SetGridValue(i, distance_grid.GridValue(i));
  }
  
  // Print statistics
  if (print_verbose) {
    printf("Updated signed distances using grid ...\n");
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

  // Return success
  return 1;
}



static int
RefineSignedDistanceUsingKdtree(R3Grid *grid, const RNArray<Point *>& points, int refinement_radius)
{
  // Start statistics
  RNTime start_time;
  start_time.Read();
  int miss_count = 0;

  // Create kdtree for point samples 
  Point tmp; int position_offset = (unsigned char *) &(tmp.position) - (unsigned char *) &tmp;
  R3Kdtree<Point *> kdtree(points, position_offset);
  
  // Compute distance to closest point sample at every negative grid cell
  RNLength grid_spacing = grid->GridToWorldScaleFactor();
  RNScalar max_distance = refinement_radius * grid_spacing;
  if (max_distance > grid->WorldBox().DiagonalLength()) max_distance = grid->WorldBox().DiagonalLength();
  if (max_distance > truncation_distance) max_distance = truncation_distance;
  for (int iz = 0; iz < grid->ZResolution(); iz++) {
    for (int iy = 0; iy < grid->YResolution(); iy++) {
      for (int ix = 0; ix < grid->XResolution(); ix++) {
        // Get current distance
        RNScalar grid_distance = grid->GridValue(ix, iy, iz);
        RNScalar sign = (grid_distance < 0) ? -1 : 1;

        // Check if in swath on negative side
        if (grid_distance >= 0) continue;
        if (grid_distance < -max_distance) continue;

        // Get world position
        R3Point world_position = grid->WorldPosition(ix, iy, iz);

        // Find closest point sample
        Point *closest = kdtree.FindClosest(world_position, 0, fabs(grid_distance) + grid_spacing);
        if (!closest) {
          // Indicate should be interpolated later
          grid->SetGridValue(ix, iy, iz, -FLT_MAX);
          miss_count++;
          continue;
        }

        // Compute distance to closest point sample
        R3Point closest_position = closest->position;
        R3Plane plane(closest->position, closest->normal);
        R3Point projected_position = world_position; projected_position.Project(plane);
        R3Vector tangent_vector = projected_position - closest->position;
        RNLength tangent_distance = tangent_vector.Length();
        if (RNIsPositive(tangent_distance)) {
          tangent_vector /= tangent_distance;
          if (tangent_distance > closest->radius) tangent_distance = closest->radius;
          closest_position = closest->position + tangent_distance * tangent_vector;
        }
        
        // Set grid value
        RNScalar closest_distance = R3Distance(world_position, closest_position);
        grid->SetGridValue(ix, iy, iz, sign * closest_distance);
      }
    }
  }

  // Fill in missing values (-FLT_MAX)
  grid->Substitute(0, 0.000123456);
  grid->Substitute(-FLT_MAX, 0);
  grid->FillHoles(grid->NEntries());
  grid->Substitute(0.000123456, 0);
  
  // Print statistics
  if (print_verbose) {
    printf("Refined signed distances using kdtree ...\n");
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
    printf("  Refinement grid radius = %d\n", refinement_radius);
    printf("  Refinement world radius = %g\n", max_distance);
    printf("  Miss count = %d\n", miss_count);
    fflush(stdout);
  }

  // Return success
  return 1;
}



////////////////////////////////////////////////////////////////////////
// Grid processing
////////////////////////////////////////////////////////////////////////

static int
SmoothDistanceAtRefinementBoundary(R3Grid *grid, int refinement_radius)
{
  // Start statistics
  RNTime start_time;
  start_time.Read();
  int count = 0;

  // Copy original grid values
  R3Grid copy_grid(*grid);
  
  // Smooth distance at grid cells near boundary between kdtree and grid estimation
  RNLength grid_spacing = grid->GridToWorldScaleFactor();
  RNScalar min_distance = (refinement_radius-2) * grid_spacing;
  RNScalar max_distance = (refinement_radius+2) * grid_spacing;
  if (min_distance < 2*grid_spacing) min_distance = 2*grid_spacing;
  for (int iz = 0; iz < grid->ZResolution(); iz++) {
    for (int iy = 0; iy < grid->YResolution(); iy++) {
      for (int ix = 0; ix < grid->XResolution(); ix++) {
        // Check if on refinement boundary
        RNScalar grid_value = copy_grid.GridValue(ix, iy, iz);
        if (fabs(grid_value) > max_distance) continue;
        if (fabs(grid_value) < min_distance) continue;

        // Compute weighted sum of distances for neighbor cells
        RNScalar sum_weight = 0;
        RNScalar sum_distance = 0;
        for (int dz = -1; dz <= 1; dz++) {
          int gz = iz + dz;
          if ((gz < 0) || (gz >= grid->ZResolution())) continue;
          for (int dy = -1; dy <= 1; dy++) {
            int gy = iy + dy;
            if ((gy < 0) || (gy >= grid->YResolution())) continue;
            for (int dx = -1; dx <= 1; dx++) {
              int gx = ix + dx;
              if ((gx < 0) || (gx >= grid->XResolution())) continue;
              RNScalar weight = pow(2, -(dx + dy + dz));
              RNScalar distance = copy_grid.GridValue(gx, gy, gz);
              sum_distance += weight * distance;
              sum_weight += weight;
            }
          }
        }
        
        // Set smoothed distance
        if (RNIsZero(sum_weight)) continue;
        RNScalar distance = sum_distance / sum_weight;
        grid->SetGridValue(ix, iy, iz, distance);
        count++;
      }
    }
  }

  // Print statistics
  if (print_verbose) {
    printf("Smoothed refinement boundary ...\n");
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
    printf("  Count = %d\n", count);
    fflush(stdout);
  }

  // Return success
  return 1;
}



static int
TruncateDistance(R3Grid *grid, RNLength truncation_distance)
{
  // Check truncation distance
  if (truncation_distance >= RN_INFINITY) return 1;

  // Start statistics
  RNTime start_time;
  start_time.Read();

  // Truncate distances
  grid->Threshold(-truncation_distance, -truncation_distance, R3_GRID_KEEP_VALUE);
  grid->Threshold(truncation_distance, R3_GRID_KEEP_VALUE, truncation_distance);

  // Print statistics
  if (print_verbose) {
    printf("Truncated distances ...\n");
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
    printf("  Threshold = %g\n", truncation_distance);
    fflush(stdout);
  }

  // Return success
  return 1;
}



////////////////////////////////////////////////////////////////////////
// Mesh extraction
////////////////////////////////////////////////////////////////////////

static RNBoolean
IsManifold(R3Mesh *mesh)
{
  // Check for boundary edges
  for (int i = 0; i < mesh->NEdges(); i++) {
    R3MeshEdge *edge = mesh->Edge(i);
    if (mesh->IsEdgeOnBoundary(edge)) return FALSE;
  }

  // Check connected components
  int ncomponents = mesh->ConnectedComponents();
  if (ncomponents != 1) return FALSE;

  // Passed all tests
  return TRUE;
}



static R3Mesh *
ExtractIsosurface(R3Grid *grid)
{
  // Start statistics
  RNTime start_time;
  start_time.Read();

  // Allocate isosurface mesh
  R3Mesh *isosurface_mesh = new R3Mesh();
  if (!isosurface_mesh) {
    RNFail("Unable to allocate isosurface mesh\n");
    return NULL;
  }

  // Extract isosurface
  // One grid_spacing outside, to account for razor-thin surfaces
  RNLength grid_spacing = grid->GridToWorldScaleFactor();
  grid->GenerateIsoSurface(grid_spacing, isosurface_mesh);

  // Print statistics
  if (print_verbose) {
    printf("Extracted isosurface ...\n");
    printf("  Time = %.2f seconds\n", start_time.Elapsed());
    printf("  # Faces = %d\n", isosurface_mesh->NFaces());
    printf("  # Edges = %d\n", isosurface_mesh->NEdges());
    printf("  # Vertices = %d\n", isosurface_mesh->NVertices());
    printf("  Manifold = %s\n", (IsManifold(isosurface_mesh)) ? "yes" : "no");
    fflush(stdout);
  }

  // Return isosurface mesh
  return isosurface_mesh;
}



////////////////////////////////////////////////////////////////////////
// Program argument parsing
////////////////////////////////////////////////////////////////////////

static int 
ParseArgs(int argc, char **argv)
{
  // Parse arguments
  argc--; argv++;
  while (argc > 0) {
    if ((*argv)[0] == '-') {
      if (!strcmp(*argv, "-v")) print_verbose = 1; 
      else if (!strcmp(*argv, "-debug")) print_debug = 1; 
      else if (!strcmp(*argv, "-estimate_sign")) estimate_sign = 1; 
      else if (!strcmp(*argv, "-estimate_sign_using_normals")) { estimate_sign = 1; sign_estimation_method = 0; }
      else if (!strcmp(*argv, "-estimate_sign_using_flood_fill")) { estimate_sign = 1; sign_estimation_method = 1; }
      else if (!strcmp(*argv, "-sign_estimation_method")) { argc--; argv++; sign_estimation_method = atoi(*argv); }
      else if (!strcmp(*argv, "-truncation_distance")) { argc--; argv++; truncation_distance = atof(*argv); }
      else if (!strcmp(*argv, "-spacing")) { argc--; argv++; target_grid_spacing = atof(*argv); }
      else if (!strcmp(*argv, "-border")) { argc--; argv++; grid_border = atof(*argv); }
      else if (!strcmp(*argv, "-max_resolution")) { argc--; argv++; grid_max_resolution = atoi(*argv); }
      else if (!strcmp(*argv, "-refinement_radius")) { argc--; argv++; refinement_radius = atoi(*argv); }
      else if (!strcmp(*argv, "-output_mesh")) { argc--; argv++; output_mesh_filename = *argv; }
      else if (!strcmp(*argv, "-output_points")) { argc--; argv++; output_points_filename = *argv; }
      else if (!strcmp(*argv, "-input_is_manifold")) { input_is_manifold = 1; }
      else if (!strcmp(*argv, "-input_is_range_scan")) { input_is_range_scan = 1; }
      else if (!strcmp(*argv, "-scan_viewpoint")) {
        argc--; argv++; scan_viewpoint[0] = atof(*argv);
        argc--; argv++; scan_viewpoint[1] = atof(*argv);
        argc--; argv++; scan_viewpoint[2] = atof(*argv);
      }
      else if (!strcmp(*argv, "-bbox")) {
        argc--; argv++; grid_bbox[0][0] = atof(*argv);
        argc--; argv++; grid_bbox[0][1] = atof(*argv);
        argc--; argv++; grid_bbox[0][2] = atof(*argv);
        argc--; argv++; grid_bbox[1][0] = atof(*argv);
        argc--; argv++; grid_bbox[1][1] = atof(*argv);
        argc--; argv++; grid_bbox[1][2] = atof(*argv);
      }
      else {
        RNFail("Invalid program argument: %s", *argv);
        exit(1);
      }
    }
    else {
      if (!input_mesh_filename) input_mesh_filename = *argv;
      else if (!output_grid_filename) output_grid_filename = *argv;
      else { RNFail("Invalid program argument: %s", *argv); exit(1); }
    }
    argv++; argc--;
  }

  // Check filenames
  if (!input_mesh_filename || !output_grid_filename) {
    RNFail("Usage: msh2df inputmeshfile outputgridfile [options]\n");
    return 0;
  }

  // Flood fill sign estimation does not work with grid distances near boundary
  if (estimate_sign && (sign_estimation_method == 1) && (refinement_radius < 3)) {
    refinement_radius = 3;
  }

  // Return OK status 
  return 1;
}



////////////////////////////////////////////////////////////////////////
// Main program
////////////////////////////////////////////////////////////////////////

int 
main(int argc, char **argv)
{
  // Parse program arguments
  if (!ParseArgs(argc, argv)) exit(-1);

  // Read mesh file
  R3Mesh *mesh = ReadMesh(input_mesh_filename);
  if (!mesh) exit(-1);

  // Update program arguments based on mesh properties
  if (grid_bbox.IsEmpty()) grid_bbox = mesh->BBox();
  if (IsManifold(mesh)) input_is_manifold = TRUE;
  if (input_is_manifold) sign_estimation_method = 1;

  // Initialize grid with appropriate dimensions
  R3Grid grid(grid_bbox, target_grid_spacing, 5, grid_max_resolution, grid_border);

  // Estimate distance
  if (!EstimateDistanceUsingGrid(&grid, mesh)) return 0;

  // Refine distance
  if (refinement_radius > 0) {
    if (!RefineDistanceUsingKdtree(&grid, mesh, refinement_radius)) return 0;
  }

  // Estimate sign
  if (estimate_sign) {
    if (!EstimateSign(&grid, mesh)) exit(-1);
  }
  
  // Update signed distance
  RNBoolean update_distances = (estimate_sign && !input_is_manifold) ? TRUE : FALSE;
  if (update_distances) {
    if (!UpdateSignedDistanceUsingGrid(&grid, mesh)) return 0;
  }

  // Create point samples
  RNArray<Point *> points;
  if (output_points_filename || update_distances) {
    if (!CreatePoints(&grid, mesh, points)) return 0;
  }

  // Write point sample file
  if (output_points_filename && !points.IsEmpty()) {
    if (!WritePoints(points, output_points_filename)) return 0;
  }
  
  // Refine signed distance
  if (update_distances && !points.IsEmpty() && (refinement_radius > 0)) {
    if (!RefineSignedDistanceUsingKdtree(&grid, points, refinement_radius)) return 0;
  }

  // Delete point samples
  if (!points.IsEmpty()) {
    if (!DeletePoints(points)) return 0;
  }

  // Truncate distance
  if (refinement_radius > 0) {
    if (!SmoothDistanceAtRefinementBoundary(&grid, refinement_radius)) exit(-1);
  }

  // Truncate distance
  if (!TruncateDistance(&grid, truncation_distance)) exit(-1);

  // Write grid
  if (output_grid_filename) {
    if (!WriteGrid(&grid, output_grid_filename)) exit(-1);
  }

  // Write mesh file
  if (output_mesh_filename) {
    // Create isosurface
    R3Mesh *isosurface = ExtractIsosurface(&grid);
    if (!isosurface) exit(-1);

    // Write mesh
    if (!WriteMesh(isosurface, output_mesh_filename)) exit(-1);

    // Delete mesh
    delete isosurface;
  }

  // Return success
  return 0;
}
