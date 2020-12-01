////////////////////////////////////////////////////////////////////////
// Source file for RGBD segmentation functions
////////////////////////////////////////////////////////////////////////



////////////////////////////////////////////////////////////////////////
// Include files
////////////////////////////////////////////////////////////////////////

#include "RGBD/RGBD.h"
#include "segmentation.h"



//////////////////////////////////////////////////////////////////////
// Namespace
////////////////////////////////////////////////////////////////////////

namespace gaps {



////////////////////////////////////////////////////////////////////////
// Parameters
////////////////////////////////////////////////////////////////////////

double max_depth = 0;
double max_neighbor_distance_factor = 16;
double max_neighbor_normal_angle = 0;
double max_neighbor_color_difference = 0;



////////////////////////////////////////////////////////////////////////
// Segmentation functions
////////////////////////////////////////////////////////////////////////

static int 
RGBDCreateSegmentationPoints(Segmentation *segmentation,
  const R2Grid& px_image, const R2Grid& py_image, const R2Grid& pz_image, 
  const R2Grid& nx_image, const R2Grid& ny_image, const R2Grid& nz_image,
  const R2Grid& depth_image, const R2Grid& radius_image,
  const R2Grid& boundary_image, const R2Image& color_image)
{
  // Allocate points 
  segmentation->point_buffer = new Point [ depth_image.NEntries() ];
  if (!segmentation->point_buffer) {
    RNFail("Unable to allocate points\n");
    return 0;
  }

  // Fill points
  for (int ix = 0; ix < depth_image.XResolution(); ix++) {
    for (int iy = 0; iy < depth_image.YResolution(); iy++) {
      int i;
      depth_image.IndicesToIndex(ix, iy, i);
      Point *point = &segmentation->point_buffer[i];

      // Check depth
      RNScalar depth = depth_image.GridValue(i);
      if (RNIsNegativeOrZero(depth)) continue;
      if ((max_depth > 0) && (depth > max_depth)) continue;
      point->depth = depth;

      // Get position
      RNScalar px = px_image.GridValue(i);
      RNScalar py = py_image.GridValue(i);
      RNScalar pz = pz_image.GridValue(i);
      point->position.Reset(px, py, pz);

      // Get normal
      RNScalar nx = nx_image.GridValue(i);
      RNScalar ny = ny_image.GridValue(i);
      RNScalar nz = nz_image.GridValue(i);
      point->normal.Reset(nx, ny, nz);

      // Get radius
      RNScalar radius = radius_image.GridValue(i);
      point->radius1 = radius;
      point->radius2 = radius;
      point->area = RN_PI * radius * radius;

      // Get color
      point->color = color_image.PixelRGB(ix, iy);
    
      // Get flags
      point->boundary = (unsigned int) (boundary_image.GridValue(i) + 0.5);

      // Set grid index
      point->data_index = i;

      // Insert point
      segmentation->points.Insert(point);
    }
  }

  // Create kdtree of points
  Point tmp; int position_offset = (unsigned char *) &(tmp.position) - (unsigned char *) &tmp;
  segmentation->kdtree = new R3Kdtree<Point *>(segmentation->points, position_offset);
  if (!segmentation->kdtree) {
    RNFail("Unable to create kdtree\n");
    return 0;
  }
  
  // Create arrays of neighbor points
  for (int i = 0; i < segmentation->points.NEntries(); i++) {
    Point *point = segmentation->points.Kth(i);
    int ix, iy, neighbor_index;
    depth_image.IndexToIndices(point->data_index, ix, iy);
    for (int s = -1; s <= 1; s++) {
      if ((ix+s < 0) || (ix+s >= depth_image.XResolution())) continue;
      for (int t = -1; t <= 1; t++) {
        if ((s == 0) && (t == 0)) continue;
        if ((iy+t < 0) || (iy+t >= depth_image.YResolution())) continue;
        depth_image.IndicesToIndex(ix+s, iy+t, neighbor_index);
        Point *neighbor = &segmentation->point_buffer[neighbor_index];

        // Check if across boundary
        if ((point->boundary & RGBD_SHADOW_BOUNDARY) && (neighbor->boundary & RGBD_SILHOUETTE_BOUNDARY)) continue;
        if ((point->boundary & RGBD_SILHOUETTE_BOUNDARY) && (neighbor->boundary & RGBD_SHADOW_BOUNDARY)) continue;

        // Check if too far away
        if (max_neighbor_distance_factor > 0) {
          RNScalar radius = (point->radius1 > neighbor->radius1) ? neighbor->radius1 : point->radius1;
          RNScalar max_neighbor_distance = (radius > 0) ? max_neighbor_distance_factor * radius : 0.25;
          RNScalar dd = R3SquaredDistance(point->position, neighbor->position);
          if (dd > max_neighbor_distance * max_neighbor_distance) continue;
        }

        // Check if too much color difference
        if (max_neighbor_color_difference > 0) {
          RNScalar dr = fabs(point->color.R() - neighbor->color.R());
          if (dr > max_neighbor_color_difference) continue;
          RNScalar dg = fabs(point->color.G() - neighbor->color.G());
          if (dg > max_neighbor_color_difference) continue;
          RNScalar db = fabs(point->color.B() - neighbor->color.B());
          if (db > max_neighbor_color_difference) continue;
        }

        // Check if too much normal angle
        if (max_neighbor_normal_angle > 0) {
          RNScalar dot = fabs(point->normal.Dot(neighbor->normal));
          RNAngle normal_angle = (dot < 1) ? RN_PI_OVER_TWO - acos(dot) : RN_PI_OVER_TWO;
          if (normal_angle > max_neighbor_normal_angle) continue;
        }
        
        // Insert neighbor
        point->neighbors.Insert(neighbor);
      }
    }
  }

  // Return success
  return 1;
}
               


static Segmentation *
RGBDCreateSegmentation(const R2Grid& px_image, const R2Grid& py_image, const R2Grid& pz_image,
  const R2Grid& nx_image, const R2Grid& ny_image, const R2Grid& nz_image,
  const R2Grid& depth_image, const R2Grid& radius_image, 
  const R2Grid& boundary_image, const R2Image& color_image,
  const R3Point& viewpoint, const R3Vector& towards, const R3Vector& up)
{
  // Adjust segmentation parameters ???
  min_cluster_points = 10 * depth_image.NEntries() / (640 * 480);
  
  // Allocate segmentation
  Segmentation *segmentation = new Segmentation();
  if (!segmentation) {
    RNFail("Unable to allocate segmentation.\n");
    return NULL;
  }

  // Create points
  if (!RGBDCreateSegmentationPoints(segmentation,
    px_image, py_image, pz_image, nx_image, ny_image, nz_image,
    depth_image, radius_image, boundary_image, color_image)) {
    RNFail("Unable to create points for segmentation.\n");
    delete segmentation;
    return 0;
  }

  // Check points
  if (segmentation->points.NEntries() == 0) {
    RNFail("Zero points for segmentation.\n");
    delete segmentation;
    return 0;
  }

  // Create clusters
  if (!segmentation->CreateClusters(PLANE_PRIMITIVE_TYPE)) {
    RNFail("Unable to create clusters for segmentation.\n");
    delete segmentation;
    return 0;
  }

  // Return segmentation
  return segmentation;
}



#if 0
static int
RGBDWriteSegmentation(Segmentation *segmentation, const char *filename)
{
  // Open file
  FILE *fp = fopen(filename, "w");
  if (!fp) {
    RNFail("Unable to open segmentation file %s\n", filename);
    return 0;
  }

  // Write clusters to file
  for (int i = 0; i < segmentation->clusters.NEntries(); i++) {
    Cluster *cluster = segmentation->clusters.Kth(i);
    fprintf(fp, "%d %d %g %g %d  %g %g %g  %g %g %g %g  %g %g %g\n",
            i+1, cluster->points.NEntries(),
            cluster->total_affinity, cluster->possible_affinity, cluster->primitive.primitive_type,
            cluster->primitive.centroid.X(), cluster->primitive.centroid.Y(), cluster->primitive.centroid.Z(),
            cluster->primitive.plane.A(), cluster->primitive.plane.B(), cluster->primitive.plane.C(), cluster->primitive.plane.D(),
            cluster->color.R(), cluster->color.G(), cluster->color.B());
  }

  // Close file
  fclose(fp);

  // Return success
  return 1;
}
#endif



int
RGBDCreateSegmentationChannel(const R2Grid& depth_image, 
  const R2Grid& px_image, const R2Grid& py_image, const R2Grid& pz_image, const R2Grid& boundary_image, 
  const R2Grid& nx_image, const R2Grid& ny_image, const R2Grid& nz_image, const R2Grid& radius_image, 
  const R2Image& color_image, R2Grid& output_segmentation_image,
  const R3Point& viewpoint, const R3Vector& towards, const R3Vector& up,
  const char *output_segmentation_filename = NULL)
{
  // Initialize segmentation image
  output_segmentation_image = depth_image;
  output_segmentation_image.Clear(R2_GRID_UNKNOWN_VALUE);
  
  // Create segmentation
  Segmentation *segmentation = RGBDCreateSegmentation(
    px_image, py_image, pz_image, nx_image, ny_image, nz_image, 
    depth_image, radius_image, boundary_image, color_image, 
    viewpoint, towards, up);

  // Check segmentation
  if (!segmentation) return 0;

  // Fill segmentation image
  for (int i = 0; i < segmentation->clusters.NEntries(); i++) {
    Cluster *cluster = segmentation->clusters.Kth(i);
    for (int j = 0; j < cluster->points.NEntries(); j++) {
      Point *point = cluster->points.Kth(j);
      if (RNIsNegativeOrZero(point->depth)) continue;
      if (point->data_index < 0) continue;
      if (point->data_index >= output_segmentation_image.NEntries()) continue;
      output_segmentation_image.SetGridValue(point->data_index, i+1);
    }
  }

  // Write segmentation file (this is a hack)
  if (output_segmentation_filename) {
    // if (!RGBDWriteSegmentation(segmentation, output_segmentation_filename)) return 0;
    if (!segmentation->WriteFile(output_segmentation_filename)) return 0;
  }
  
  // Delete segmentation
  delete segmentation;

  // Return success
  return 1;
}




} // end namespace
