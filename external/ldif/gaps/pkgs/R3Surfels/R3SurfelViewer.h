/* Include file for the R3 surfel viewer class */
#ifndef __R3__SURFEL__VIEWER__H__
#define __R3__SURFEL__VIEWER__H__



////////////////////////////////////////////////////////////////////////
// NAMESPACE 
////////////////////////////////////////////////////////////////////////

namespace gaps {



////////////////////////////////////////////////////////////////////////
// CLASS DEFINITION
////////////////////////////////////////////////////////////////////////

class R3SurfelViewer {
public:
  //////////////////////////////////////////
  //// CONSTRUCTOR/DESTRUCTOR FUNCTIONS ////
  //////////////////////////////////////////

  // Constructor function
  R3SurfelViewer(R3SurfelScene *scene);

  // Destructor function
  virtual ~R3SurfelViewer(void);


  ////////////////////////////////////
  //// UI EVENT HANDLER FUNCTIONS ////
  ////////////////////////////////////

  // Call these at beginning and end of execution
  virtual void Initialize(void);
  virtual void Terminate(void);

  // Call these to handle user input events
  virtual int Redraw(void);
  virtual int Resize(int width, int height);
  virtual int MouseMotion(int x, int y);
  virtual int MouseButton(int x, int y, int button, int state, int shift, int ctrl, int alt);
  virtual int Keyboard(int x, int y, int key, int shift, int ctrl, int alt);
  virtual int Idle(void);


  ////////////////////////////////
  //// SCENE ACCESS FUNCTIONS ////
  ////////////////////////////////

  // Scene properties
  R3SurfelScene *Scene(void) const;


  ///////////////////////////////////
  //// PROPERTY ACCESS FUNCTIONS ////
  ///////////////////////////////////

  // Viewing properties
  const R3Camera& Camera(void) const;
  const R2Viewport& Viewport(void) const;
  const R3Box& ViewingExtent(void) const;
  const R3Point& CenterPoint(void) const;
  RNLength ImagePlaneDepth(void) const;
  RNScalar SurfelSize(void) const;

  // Visibility properties (0=off, 1=on)
  int SurfelVisibility(void) const;
  int NormalVisibility(void) const;
  int BackfacingVisibility(void) const;
  int AerialVisibility(void) const;
  int TerrestrialVisibility(void) const;
  int HumanLabeledObjectVisibility(void) const;
  int ObjectPropertyVisibility(void) const;
  int ObjectRelationshipVisibility(void) const;
  int NodeBBoxVisibility(void) const;
  int BlockBBoxVisibility(void) const;
  int ScanViewpointVisibility(void) const;
  int ImageViewpointVisibility(void) const;
  int ImagePlaneVisibility(void) const;
  int ImageInsetVisibility(void) const;
  int ImagePointsVisibility(void) const;
  int CenterPointVisibility(void) const;
  int ViewingExtentVisibility(void) const;
  int AxesVisibility(void) const;

  // Color properties
  int SurfelColorScheme(void) const;
  const char *SurfelColorSchemeName(void) const;
  const RNRgb& NormalColor(void) const;
  const RNRgb& BackgroundColor(void) const;
  const RNRgb& ObjectPropertyColor(void) const;
  const RNRgb& NodeBBoxColor(void) const;
  const RNRgb& BlockBBoxColor(void) const;
  const RNRgb& ScanViewpointColor(void) const;
  const RNRgb& ImageViewpointColor(void) const;
  const RNRgb& CenterPointColor(void) const;

  // Viewing parameters
  RNScalar TargetResolution(void) const;
  RNScalar FocusRadius(void) const;

  // Selection properties
  R3SurfelPoint *SelectedPoint(void) const;
  R3SurfelImage *SelectedImage(void) const;
  
  // Frame rate statistics
  RNScalar FrameRate(void) const;
  RNScalar FrameTime(void) const;
  RNScalar CurrentTime(void) const;

  // Image loading properties
  const char *ImageDirectory(void) const;
  
  
  /////////////////////////////////////////
  //// PROPERTY MANIPULATION FUNCTIONS ////
  /////////////////////////////////////////

  // Viewing property manipulation
  void ResetCamera(void);
  void SetCamera(const R3Camera& camera);
  void ZoomCamera(RNScalar scale = 10);
  void SetViewport(const R2Viewport& viewport);
  void SetViewingExtent(const R3Box& box);
  void SetCenterPoint(const R3Point& point);
  void SetImagePlaneDepth(RNLength depth);
  void SetSurfelSize(RNScalar npixels);

  // Visibility manipulation (0=off, 1=on, -1=toggle)
  void SetSurfelVisibility(int visibility);
  void SetNormalVisibility(int visibility);
  void SetBackfacingVisibility(int visibility);
  void SetAerialVisibility(int visibility);
  void SetTerrestrialVisibility(int visibility);
  void SetHumanLabeledObjectVisibility(int visibility);
  void SetObjectPropertyVisibility(int visibility);
  void SetObjectRelationshipVisibility(int visibility);
  void SetNodeBBoxVisibility(int visibility);
  void SetBlockBBoxVisibility(int visibility);
  void SetScanViewpointVisibility(int visibility);
  void SetImageViewpointVisibility(int visibility);
  void SetImagePlaneVisibility(int visibility);
  void SetImageInsetVisibility(int visibility);
  void SetImagePointsVisibility(int visibility);
  void SetCenterPointVisibility(int visibility);
  void SetViewingExtentVisibility(int visibility);
  void SetAxesVisibility(int visibility);

  // Color manipulation
  void SetSurfelColorScheme(int scheme);
  void SetNormalColor(const RNRgb& color);
  void SetBackgroundColor(const RNRgb& color);
  void SetObjectPropertyColor(const RNRgb& color);
  void SetNodeBBoxColor(const RNRgb& color);
  void SetBlockBBoxColor(const RNRgb& color);
  void SetScanViewpointColor(const RNRgb& color);
  void SetImageViewpointColor(const RNRgb& color);
  void SetCenterPointColor(const RNRgb& color);

  // Selection manipulation
  void SelectPoint(R3SurfelBlock *block, int surfel_index);
  void SelectImage(R3SurfelImage *image,
    RNBoolean update_working_set = TRUE,
    RNBoolean jump_to_viewpoint = FALSE);

  // Working set parameters
  void SetTargetResolution(RNScalar resolution);
  void SetFocusRadius(RNScalar radius);

  // Image loading
  void SetImageDirectory(const char *directory_name);

  // Image input/output
  int WriteImage(const char *filename);


  /////////////////////////////////////////
  //// PICKING FUNCTIONS               ////
  /////////////////////////////////////////

  // Pick utility functions
  R3SurfelImage *PickImage(int x, int y, R3Point *picked_position = NULL);
  R3SurfelNode *PickNode(int x, int y, R3Point *picked_position = NULL,
    R3SurfelBlock **picked_block = NULL, int *picked_surfel_index = NULL,
    RNBoolean exclude_nonobjects = FALSE);


  /////////////////////////////////////////
  //// DRAWING UTILITY FUNCTIONS       ////
  /////////////////////////////////////////

  void LoadColor(int k) const;
  void LoadColor(double value) const;
  void EnableViewingExtent(void) const;
  void DisableViewingExtent(void) const;
  void DrawViewingExtent(void) const;
  void DrawNode(R3SurfelNode *node, RNFlags color_draw_flags = 0) const;
  int CheckNodeVisibility(R3SurfelNode *node) const;
  
  
////////////////////////////////////////////////////////////////////////
// INTERNAL STUFF BELOW HERE
////////////////////////////////////////////////////////////////////////

  // Object editing 
  int SplitLeafNodes(R3SurfelNode *source_node, const R3SurfelConstraint& constraint, 
    RNArray<R3SurfelNode *> *nodesA = NULL, RNArray<R3SurfelNode *> *nodesB = NULL);
  int SplitObject(R3SurfelObject *object, const R3SurfelConstraint& constraint,
    R3SurfelObject **objectA = NULL, R3SurfelObject **objectB = NULL);

  // Working set management
  virtual void EmptyWorkingSet(void);
  virtual void UpdateWorkingSet(void);
  virtual void UpdateWorkingSet(const R3Viewer& view);
  virtual void UpdateWorkingSet(const R3Point& center, RNScalar target_resolution, RNScalar focus_radius);
  virtual void InsertIntoWorkingSet(R3SurfelNode *node, RNBoolean full_resolution = FALSE);
  virtual void RemoveFromWorkingSet(R3SurfelNode *node, RNBoolean full_resolution = FALSE);

protected:
  // Tree manipulation
  void SetScene(R3SurfelScene *scene);

  // Viewing utility functions
  void RotateWorld(RNScalar factor, const R3Point& origin, int, int, int dx, int dy);

  // Draw functions
  void DrawObject(R3SurfelObject *object, RNFlags flags = R3_SURFEL_DEFAULT_DRAW_FLAGS) const;

  // Management functions
  void ReadCoarsestBlocks(RNScalar max_complexity);
  void ReleaseCoarsestBlocks(RNScalar max_complexity);


protected:
  // Tree properties
  R3SurfelScene *scene;

  // Node working set
  R3SurfelNodeSet resident_nodes;

  // Viewing properties
  R3Viewer viewer;
  R3Box viewing_extent;
  R3Point center_point;
  R3SurfelPoint *selected_point;
  R3SurfelImage *selected_image;
  R2Texture current_image_texture;
  RNLength image_plane_depth;
  RNScalar surfel_size;

  // Visibility properties
  int surfel_visibility;
  int normal_visibility;
  int backfacing_visibility;
  int aerial_visibility;
  int terrestrial_visibility;
  int human_labeled_object_visibility;
  int object_property_visibility;
  int object_relationship_visibility;
  int node_bbox_visibility;
  int block_bbox_visibility;
  int scan_viewpoint_visibility;
  int image_viewpoint_visibility;
  int image_plane_visibility;
  int image_inset_visibility;
  int image_points_visibility;
  int center_point_visibility;
  int viewing_extent_visibility;
  int axes_visibility;

  // Color properties
  int surfel_color_scheme;
  RNRgb normal_color;
  RNRgb background_color;
  RNRgb object_property_color;
  RNRgb node_bbox_color;
  RNRgb block_bbox_color;
  RNRgb scan_viewpoint_color;
  RNRgb image_viewpoint_color;
  RNRgb center_point_color;

  // Shape primitive properties
  RNFlags shape_draw_flags;
  
  // Working set parameters
  RNBoolean adapt_working_set_automatically;
  RNScalar target_resolution;
  RNScalar focus_radius;

  // UI state
  int window_height;
  int window_width;
  int mouse_button[3];
  int mouse_position[2];
  int mouse_down_position[2];
  int mouse_drag_distance_squared;
  int shift_down;
  int ctrl_down;
  int alt_down;

  // Timing
  RNTime start_timer;
  RNTime frame_timer;
  RNScalar frame_time;

  // Screen capture
  char *screenshot_name;
};



////////////////////////////////////////////////////////////////////////
// KEY CONSTANT DEFINITIONS
////////////////////////////////////////////////////////////////////////

#define R3_SURFEL_VIEWER_ESC_KEY           27
#define R3_SURFEL_VIEWER_SPACE_KEY         32
#define R3_SURFEL_VIEWER_DEL_KEY          127

#define R3_SURFEL_VIEWER_LEFT_KEY        1024
#define R3_SURFEL_VIEWER_RIGHT_KEY       1025
#define R3_SURFEL_VIEWER_DOWN_KEY        1026
#define R3_SURFEL_VIEWER_UP_KEY          1027
#define R3_SURFEL_VIEWER_HOME_KEY        1028
#define R3_SURFEL_VIEWER_END_KEY         1029
#define R3_SURFEL_VIEWER_INSERT_KEY      1030
#define R3_SURFEL_VIEWER_PAGE_DOWN_KEY   1031
#define R3_SURFEL_VIEWER_PAGE_UP_KEY     1032

#define R3_SURFEL_VIEWER_F1_KEY          1064
#define R3_SURFEL_VIEWER_F2_KEY          1065
#define R3_SURFEL_VIEWER_F3_KEY          1066
#define R3_SURFEL_VIEWER_F4_KEY          1067
#define R3_SURFEL_VIEWER_F5_KEY          1068
#define R3_SURFEL_VIEWER_F6_KEY          1069
#define R3_SURFEL_VIEWER_F7_KEY          1070
#define R3_SURFEL_VIEWER_F8_KEY          1071
#define R3_SURFEL_VIEWER_F9_KEY          1072
#define R3_SURFEL_VIEWER_F10_KEY         1073
#define R3_SURFEL_VIEWER_F11_KEY         1074
#define R3_SURFEL_VIEWER_F12_KEY         1075



////////////////////////////////////////////////////////////////////////
// COLOR SCHEME CONSTANT DEFINITIONS
////////////////////////////////////////////////////////////////////////

enum {
  R3_SURFEL_VIEWER_COLOR_BY_RGB,
  R3_SURFEL_VIEWER_COLOR_BY_SHADING,
  R3_SURFEL_VIEWER_COLOR_BY_HEIGHT,
  R3_SURFEL_VIEWER_COLOR_BY_NORMAL,
  R3_SURFEL_VIEWER_COLOR_BY_SCAN,
  R3_SURFEL_VIEWER_COLOR_BY_OBJECT,
  R3_SURFEL_VIEWER_COLOR_BY_NODE,
  R3_SURFEL_VIEWER_COLOR_BY_BLOCK,
  R3_SURFEL_VIEWER_COLOR_BY_CURRENT_LABEL,
  R3_SURFEL_VIEWER_COLOR_BY_GROUND_TRUTH_LABEL,
  R3_SURFEL_VIEWER_COLOR_BY_VIEWPOINT,
  R3_SURFEL_VIEWER_NUM_COLOR_SCHEMES
};



////////////////////////////////////////////////////////////////////////
// INLINE FUNCTIONS
////////////////////////////////////////////////////////////////////////

inline const R3Camera& R3SurfelViewer::
Camera(void) const
{
  // Return camera
  return viewer.Camera();
}



inline const R2Viewport& R3SurfelViewer::
Viewport(void) const
{
  // Return viewport
  return viewer.Viewport();
}



inline const R3Box& R3SurfelViewer::
ViewingExtent(void) const
{
  // Return viewing extent
  return viewing_extent;
}



inline const R3Point& R3SurfelViewer::
CenterPoint(void) const
{
  // Return center point
  return center_point;
}



inline RNLength R3SurfelViewer::
ImagePlaneDepth(void) const
{
  // Return depth at which image plane will be drawn
  return image_plane_depth;
}



inline RNScalar R3SurfelViewer::
SurfelSize(void) const
{
  // Return surfel size
  return surfel_size;
}



inline int R3SurfelViewer::
SurfelVisibility(void) const
{
  // Return surfel visibililty
  return surfel_visibility;
}



inline int R3SurfelViewer::
NormalVisibility(void) const
{
  // Return normal visibililty
  return normal_visibility;
}



inline int R3SurfelViewer::
BackfacingVisibility(void) const
{
  // Return backfacing visibililty
  return backfacing_visibility;
}



inline int R3SurfelViewer::
AerialVisibility(void) const
{
  // Return background visibililty
  return aerial_visibility;
}



inline int R3SurfelViewer::
TerrestrialVisibility(void) const
{
  // Return background visibililty
  return terrestrial_visibility;
}



inline int R3SurfelViewer::
HumanLabeledObjectVisibility(void) const
{
  // Return human labeled object visibililty
  return human_labeled_object_visibility;
}



inline int R3SurfelViewer::
ObjectPropertyVisibility(void) const
{
  // Return object property visibililty
  return object_property_visibility;
}



inline int R3SurfelViewer::
ObjectRelationshipVisibility(void) const
{
  // Return object relationship visibililty
  return object_relationship_visibility;
}



inline int R3SurfelViewer::
NodeBBoxVisibility(void) const
{
  // Return node bbox visibililty
  return node_bbox_visibility;
}



inline int R3SurfelViewer::
BlockBBoxVisibility(void) const
{
  // Return block bbox visibililty
  return block_bbox_visibility;
}



inline int R3SurfelViewer::
ScanViewpointVisibility(void) const
{
  // Return scan viewpoint visibililty
  return scan_viewpoint_visibility;
}



inline int R3SurfelViewer::
ImageViewpointVisibility(void) const
{
  // Return image viewpoint visibililty
  return image_viewpoint_visibility;
}



inline int R3SurfelViewer::
ImagePlaneVisibility(void) const
{
  // Return image plane visibililty
  return image_plane_visibility;
}



inline int R3SurfelViewer::
ImageInsetVisibility(void) const
{
  // Return image inset visibililty
  return image_inset_visibility;
}



inline int R3SurfelViewer::
ImagePointsVisibility(void) const
{
  // Return image points visibililty
  return image_points_visibility;
}



inline int R3SurfelViewer::
CenterPointVisibility(void) const
{
  // Return center point visibililty
  return center_point_visibility;
}



inline int R3SurfelViewer::
AxesVisibility(void) const
{
  // Return axes visibililty
  return axes_visibility;
}



inline int R3SurfelViewer::
ViewingExtentVisibility(void) const
{
  // Return viewing extent visibililty
  return viewing_extent_visibility;
}



inline int R3SurfelViewer::
SurfelColorScheme(void) const
{
  // Return color scheme for drawing surfels
  return surfel_color_scheme;
}



inline const RNRgb& R3SurfelViewer::
NormalColor(void) const
{
  // Return normal color
  return normal_color;
}



inline const RNRgb& R3SurfelViewer::
BackgroundColor(void) const
{
  // Return background color
  return background_color;
}



inline const RNRgb& R3SurfelViewer::
ObjectPropertyColor(void) const
{
  // Return object property color
  return object_property_color;
}



inline const RNRgb& R3SurfelViewer::
NodeBBoxColor(void) const
{
  // Return node bbox color
  return node_bbox_color;
}



inline const RNRgb& R3SurfelViewer::
BlockBBoxColor(void) const
{
  // Return block bbox color
  return block_bbox_color;
}



inline const RNRgb& R3SurfelViewer::
ScanViewpointColor(void) const
{
  // Return scan viewpoint color
  return scan_viewpoint_color;
}



inline const RNRgb& R3SurfelViewer::
ImageViewpointColor(void) const
{
  // Return image viewpoint color
  return image_viewpoint_color;
}



inline const RNRgb& R3SurfelViewer::
CenterPointColor(void) const
{
  // Return center point color
  return center_point_color;
}



inline RNScalar R3SurfelViewer::
TargetResolution(void) const
{
  // Return target resolution
  return target_resolution;
}



inline RNScalar R3SurfelViewer::
FocusRadius(void) const
{
  // Return focus radius
  return focus_radius;
}



inline RNScalar R3SurfelViewer::
FrameRate(void) const
{
  // Return frame rate
  if (frame_time == 0) return 0;
  else return 1.0 / frame_time;
}



inline RNScalar R3SurfelViewer::
FrameTime(void) const
{
  // Return number of seconds per frame refresh
  return frame_time;
}



inline RNScalar R3SurfelViewer::
CurrentTime(void) const
{
  // Return number of seconds since start up
  return start_timer.Elapsed();
}



inline R3SurfelScene *R3SurfelViewer::
Scene(void) const
{
  // Return scene
  return scene;
}


inline void R3SurfelViewer::
SetCamera(const R3Camera& camera)
{
  // Set camera
  viewer.SetCamera(camera);
}



inline void R3SurfelViewer::
SetViewport(const R2Viewport& viewport)
{
  // Set viewport
  viewer.SetViewport(viewport);
}



inline void R3SurfelViewer::
SetViewingExtent(const R3Box& box)
{
  // Set viewing extent
  viewing_extent = box;
}



inline void R3SurfelViewer::
SetCenterPoint(const R3Point& point)
{
  // Set center point
  center_point = point;

  // Update working set
  UpdateWorkingSet(center_point, target_resolution, focus_radius);
}



inline void R3SurfelViewer::
SetImagePlaneDepth(RNLength depth)
{
  // Set depth at which image plane will be drawn
  this->image_plane_depth = depth;
}



inline void R3SurfelViewer::
SetSurfelSize(RNScalar npixels)
{
  // Set surfel size
  this->surfel_size = npixels;
}



inline void R3SurfelViewer::
SetSurfelVisibility(int visibility)
{
  // Set surfel visibililty
  if (visibility == -1) surfel_visibility = 1 - surfel_visibility;
  else if (visibility == 0) surfel_visibility = 0;
  else surfel_visibility = 1;
}



inline void R3SurfelViewer::
SetNormalVisibility(int visibility)
{
  // Set normal visibililty
  if (visibility == -1) normal_visibility = 1 - normal_visibility;
  else if (visibility == 0) normal_visibility = 0;
  else normal_visibility = 1;
}



inline void R3SurfelViewer::
SetBackfacingVisibility(int visibility)
{
  // Set backfacing visibililty
  if (visibility == -1) backfacing_visibility = 1 - backfacing_visibility;
  else if (visibility == 0) backfacing_visibility = 0;
  else backfacing_visibility = 1;
}



inline void R3SurfelViewer::
SetAerialVisibility(int visibility)
{
  // Set background visibililty
  if (visibility == -1) aerial_visibility = 1 - aerial_visibility;
  else if (visibility == 0) aerial_visibility = 0;
  else aerial_visibility = 1;
}



inline void R3SurfelViewer::
SetTerrestrialVisibility(int visibility)
{
  // Set background visibililty
  if (visibility == -1) terrestrial_visibility = 1 - terrestrial_visibility;
  else if (visibility == 0) terrestrial_visibility = 0;
  else terrestrial_visibility = 1;
}



inline void R3SurfelViewer::
SetHumanLabeledObjectVisibility(int visibility)
{
  // Set human labeled object visibililty
  if (visibility == -1) human_labeled_object_visibility = 1 - human_labeled_object_visibility;
  else if (visibility == 0) human_labeled_object_visibility = 0;
  else human_labeled_object_visibility = 1;
}



inline void R3SurfelViewer::
SetObjectPropertyVisibility(int visibility)
{
  // Set object property visibililty
  if (visibility == -1) object_property_visibility = 1 - object_property_visibility;
  else if (visibility == 0) object_property_visibility = 0;
  else object_property_visibility = 1;
}



inline void R3SurfelViewer::
SetObjectRelationshipVisibility(int visibility)
{
  // Set object relationship visibililty
  if (visibility == -1) object_relationship_visibility = 1 - object_relationship_visibility;
  else if (visibility == 0) object_relationship_visibility = 0;
  else object_relationship_visibility = 1;
}



inline void R3SurfelViewer::
SetNodeBBoxVisibility(int visibility)
{
  // Set node bbox visibililty
  if (visibility == -1) node_bbox_visibility = 1 - node_bbox_visibility;
  else if (visibility == 0) node_bbox_visibility = 0;
  else node_bbox_visibility = 1;
}



inline void R3SurfelViewer::
SetBlockBBoxVisibility(int visibility)
{
  // Set block bbox visibililty
  if (visibility == -1) block_bbox_visibility = 1 - block_bbox_visibility;
  else if (visibility == 0) block_bbox_visibility = 0;
  else block_bbox_visibility = 1;
}



inline void R3SurfelViewer::
SetScanViewpointVisibility(int visibility)
{
  // Set scan viewpoint visibililty
  if (visibility == -1) scan_viewpoint_visibility = 1 - scan_viewpoint_visibility;
  else if (visibility == 0) scan_viewpoint_visibility = 0;
  else scan_viewpoint_visibility = 1;
}



inline void R3SurfelViewer::
SetImageViewpointVisibility(int visibility)
{
  // Set image viewpoint visibililty
  if (visibility == -1) image_viewpoint_visibility = 1 - image_viewpoint_visibility;
  else if (visibility == 0) image_viewpoint_visibility = 0;
  else image_viewpoint_visibility = 1;
}



inline void R3SurfelViewer::
SetImagePlaneVisibility(int visibility)
{
  // Set image plane visibililty
  if (visibility == -1) image_plane_visibility = 1 - image_plane_visibility;
  else if (visibility == 0) image_plane_visibility = 0;
  else image_plane_visibility = 1;
}



inline void R3SurfelViewer::
SetImageInsetVisibility(int visibility)
{
  // Set image inset visibililty
  if (visibility == -1) image_inset_visibility = 1 - image_inset_visibility;
  else if (visibility == 0) image_inset_visibility = 0;
  else image_inset_visibility = 1;
}



inline void R3SurfelViewer::
SetImagePointsVisibility(int visibility)
{
  // Set image points visibililty
  if (visibility == -1) image_points_visibility = 1 - image_points_visibility;
  else if (visibility == 0) image_points_visibility = 0;
  else image_points_visibility = 1;
}



inline void R3SurfelViewer::
SetCenterPointVisibility(int visibility)
{
  // Set center point visibililty
  if (visibility == -1) center_point_visibility = 1 - center_point_visibility;
  else if (visibility == 0) center_point_visibility = 0;
  else center_point_visibility = 1;
}



inline void R3SurfelViewer::
SetViewingExtentVisibility(int visibility)
{
  // Set viewing extent visibililty
  if (visibility == -1) viewing_extent_visibility = 1 - viewing_extent_visibility;
  else if (visibility == 0) viewing_extent_visibility = 0;
  else viewing_extent_visibility = 1;
}



inline void R3SurfelViewer::
SetAxesVisibility(int visibility)
{
  // Set axes visibililty
  if (visibility == -1) axes_visibility = 1 - axes_visibility;
  else if (visibility == 0) axes_visibility = 0;
  else axes_visibility = 1;
}



inline void R3SurfelViewer::
SetSurfelColorScheme(int scheme)
{
  // Set color scheme for drawing surfels
  surfel_color_scheme = scheme;
}



inline void R3SurfelViewer::
SetNormalColor(const RNRgb& color)
{
  // Set normal color
  normal_color = color;
}



inline void R3SurfelViewer::
SetBackgroundColor(const RNRgb& color)
{
  // Set background color
  background_color = color;
}



inline void R3SurfelViewer::
SetObjectPropertyColor(const RNRgb& color)
{
  // Set object property color
  object_property_color = color;
}



inline void R3SurfelViewer::
SetNodeBBoxColor(const RNRgb& color)
{
  // Set node bbox color
  node_bbox_color = color;
}



inline void R3SurfelViewer::
SetBlockBBoxColor(const RNRgb& color)
{
  // Set block bbox color
  block_bbox_color = color;
}



inline void R3SurfelViewer::
SetScanViewpointColor(const RNRgb& color)
{
  // Set scan viewpoint color
  scan_viewpoint_color = color;
}



inline void R3SurfelViewer::
SetImageViewpointColor(const RNRgb& color)
{
  // Set image viewpoint color
  image_viewpoint_color = color;
}



inline void R3SurfelViewer::
SetCenterPointColor(const RNRgb& color)
{
  // Set center point color
  center_point_color = color;
}



inline void R3SurfelViewer::
SetTargetResolution(RNScalar target_resolution)
{
  // Set target resolution
  this->target_resolution = target_resolution;

  // Update working set
  UpdateWorkingSet(center_point, target_resolution, focus_radius);
}



inline void R3SurfelViewer::
SetFocusRadius(RNScalar focus_radius)
{
  // Set focus radius
  this->focus_radius = focus_radius;

  // Update working set
  UpdateWorkingSet(center_point, target_resolution, focus_radius);
}



inline R3SurfelPoint *R3SurfelViewer::
SelectedPoint(void) const
{
  // Return currently selected point
  return selected_point;
}



inline R3SurfelImage *R3SurfelViewer::
SelectedImage(void) const
{
  // Return currently selected image
  return selected_image;
}



inline void R3SurfelViewer::
UpdateWorkingSet(void)
{
  // Update working set
  UpdateWorkingSet(center_point, target_resolution, focus_radius);
}



inline void
LoadUnsignedInt(unsigned int value)
{
  // Load identifer
  unsigned char rgba[4];
  rgba[0] = value & 0xFF;
  rgba[1] = (value >> 8) & 0xFF;
  rgba[2] = (value >> 16) & 0xFF;
  rgba[3] = (value >> 24) & 0xFF;
  glColor4ubv(rgba);
}



inline unsigned int
DecodeUnsignedInt(unsigned char rgba[4])
{
  // Decode identifer
  int r = rgba[0] & 0xFF;
  int g = rgba[1] & 0xFF;
  int b = rgba[2] & 0xFF;
  int a = rgba[3] & 0xFF;
  int value = r | (g << 8) | (b << 16) | (a << 24);
  return value;
}



// End namespace
}



#endif
