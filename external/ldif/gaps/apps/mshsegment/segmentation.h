//////////////////////////////////////////////////////////////////////
// Include file for code that does segmentation of points
////////////////////////////////////////////////////////////////////////
#ifndef __R3__SEGMENTATION__H__
#define __R3__SEGMENTATION__H__



//////////////////////////////////////////////////////////////////////
// Include files
////////////////////////////////////////////////////////////////////////

#include "R3Shapes/R3Shapes.h"



//////////////////////////////////////////////////////////////////////
// Namespace
////////////////////////////////////////////////////////////////////////

namespace gaps {



//////////////////////////////////////////////////////////////////////
// Parameters
////////////////////////////////////////////////////////////////////////

extern int min_clusters;
extern int max_clusters;
extern int min_cluster_points;
extern double min_cluster_area;
extern double min_cluster_coverage;
extern double max_cluster_diameter;
extern double max_cluster_primitive_distance;
extern double max_cluster_normal_angle;
extern double max_cluster_color_difference;
extern double max_cluster_timestamp_difference;
extern double max_pair_centroid_distance;
extern double max_pair_primitive_distance;
extern double max_pair_normal_angle;
extern double max_pair_color_difference;
extern double max_pair_timestamp_difference;
extern double min_pair_affinity;
extern int max_refinement_iterations;  
extern int max_reassignment_iterations;  
extern RNBoolean equalize_cluster_sizes;
extern RNBoolean balance_cluster_sizes;
extern RNBoolean favor_convex_clusters;
extern RNBoolean initialize_hierarchically;
extern RNBoolean allow_outlier_points;
extern RNBoolean refine_boundaries;
extern int print_progress;




///////////////////////////////////////////////////////////////////////
// Shape types
////////////////////////////////////////////////////////////////////////

enum {
  NULL_PRIMITIVE_TYPE,
  POINT_PRIMITIVE_TYPE,
  LINE_PRIMITIVE_TYPE,
  PLANE_PRIMITIVE_TYPE,
  PLANAR_GRID_PRIMITIVE_TYPE,
  NUM_PRIMITIVE_TYPES
};



////////////////////////////////////////////////////////////////////////
// Type definitions
////////////////////////////////////////////////////////////////////////

struct Point {
public:
  Point(void);
public:
  RNScalar depth;
  R3Point position;
  R3Vector normal;
  R3Vector tangent;
  RNLength radius1;
  RNLength radius2;
  RNScalar timestamp;
  unsigned int identifier;
  RNArea area;
  RNRgb color;
  unsigned int boundary;
  RNArray<Point *> neighbors;
  struct Cluster *cluster;
  RNScalar cluster_affinity;
  int cluster_index;
  int data_index;
  int mark;
};

struct Primitive {
  Primitive(int primitive_type = 0);
  Primitive(const Primitive& primitive);
  Primitive(Point *seed_point, const RNArray<Point *> *points = NULL);
  RNLength Distance(const R3Point& position) const;
  void Update(const R3Point& point);
  void Update(const R3Line& line);
  void Update(const R3Plane& plane);
  void Update(Point *seed_point = NULL, const RNArray<Point *> *points = NULL);
  void Update(Primitive primitive1, Primitive primitive2, RNScalar weight1 = 1.0, RNScalar weight2 = 1.0);
public:
  int primitive_type;
  R3Box bbox;
  R3Point centroid;
  R3Line line;
  R3Plane plane;
};

struct Cluster {
public:
  Cluster(Point *seed_point = NULL, int primitive_type = 0);
  Cluster(Point *seed_point, const Primitive& primitive);
  Cluster(Cluster *child1, Cluster *child2);
  ~Cluster(void);
  RNScalar Coverage(void);
  R3Triad PrincipleAxes(R3Point *returned_center = NULL, RNScalar *returned_variances = NULL) const;
  void EmptyPoints(void);
  void InsertPoint(Point *point, RNScalar affinity = 1.0);
  void RemovePoint(Point *point);
  void InsertChild(Cluster *child);
  void RemoveChild(Cluster *child);
  int UpdatePoints(const R3Kdtree<Point *> *kdtree);
  int UpdatePrimitive(void);
  int UpdateColor(void);
  int UpdateTimestamp(void);
  int UpdateArea(void);
  RNScalar Affinity(Point *point) const;
  RNScalar Affinity(Cluster *cluster) const;
public:
  Point *seed_point;
  RNArray<Point *> points;
  Cluster *parent;
  RNArray<Cluster *> children;
  RNArray<struct Pair *> pairs;
  Primitive primitive;
  RNArea area;
  RNRgb color;
  RNScalar timestamp;
  RNScalar possible_affinity; 
  RNScalar total_affinity;
  struct Segmentation *segmentation;
  int segmentation_index;
};

struct Pair {
public:
  Pair(Cluster *cluster1 = NULL, Cluster *cluster2 = NULL, RNScalar affinity = 0);
  ~Pair(void);
public:
  Cluster *clusters[2];
  int cluster_index[2];
  RNScalar affinity; 
  Pair **heapentry;
};

struct Segmentation {
public:
  Segmentation(void);
  ~Segmentation(void);
  RNScalar Affinity(void) const;
  int NUnclusteredPoints(void) const;
public:
  int CreateNeighbors(int max_neighbor_count = 16,
    double max_neighbor_distance = 0,
    double max_neighbor_primitive_distance = 0.01,
    double max_neighbor_normal_angle = RN_PI / 4.0,
    double max_neighbor_color_difference = 0,
    double max_neighbor_distance_factor = 10,
    double max_timestamp_difference = 0);
public:
  int AssignPoints(void);  
  int CreateClusters(int primitive_type);
  int CreateSingletonClusters(int primitive_type);
  int CreateRegionGrowingClusters(int primitive_type);
  int RefineClusters(void);  
  int ReassignClusters(void);  
  int DeleteClusters(void);  
  int MergeClusters(void);  
  int MergeSmallClusters(void);  
  int SplitClusters(void);
  int RefineBoundaries(void);
  int WriteFile(const char *filename) const;
public:
  RNArray<Point *> points;
  R3Kdtree<Point *> *kdtree;
  RNArray<Cluster *> clusters;
  Point *point_buffer;
};



} // end namespace

  

#endif
