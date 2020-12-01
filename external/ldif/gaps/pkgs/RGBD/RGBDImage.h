////////////////////////////////////////////////////////////////////////
// Include file for RGBDImage class
////////////////////////////////////////////////////////////////////////

#ifndef __RGBD__IMAGE__H__
#define __RGBD__IMAGE__H__



////////////////////////////////////////////////////////////////////////
// Namespace
////////////////////////////////////////////////////////////////////////

namespace gaps {

  

////////////////////////////////////////////////////////////////////////
// Class definition
////////////////////////////////////////////////////////////////////////

class RGBDImage {
public:
  // Constructors/destructors
  RGBDImage(void);
  RGBDImage(const char *color_filename, const char *depth_filename,
    const R3Matrix& intrinsics_matrix, const R4Matrix& camera_to_world_matrix,
    int width = 0, int height = 0, RNScalar timestamp = 0);
  virtual ~RGBDImage(void);

  // Configuration access functions
  RGBDConfiguration *Configuration(void) const;
  int ConfigurationIndex(void) const;

  // Channel access functions
  int NChannels(void) const;
  R2Grid *Channel(int channel_index) const;
  R2Grid *RedChannel(void) const;
  R2Grid *GreenChannel(void) const;
  R2Grid *BlueChannel(void) const;
  R2Grid *DepthChannel(void) const;
  R2Grid *CategoryChannel(void) const;
  R2Grid *InstanceChannel(void) const;
  
  // Pixel property access functions
  int NPixels(int dim) const;
  RNRgb PixelColor(int ix, int iy) const;
  RNScalar PixelDepth(int ix, int iy) const;
  RNScalar PixelCategory(int ix, int iy) const;
  RNScalar PixelInstance(int ix, int iy) const;
  R3Point PixelWorldPosition(int ix, int iy) const;
  R3Vector PixelWorldNormal(int ix, int iy) const;
  RNScalar PixelChannelValue(int ix, int iy, int channel_index) const;

  // Point property access functions
  RNRgb PixelColor(const R2Point& image_position) const;
  RNScalar PixelDepth(const R2Point& image_position) const;
  R3Point PixelWorldPosition(const R2Point& image_position) const;
  R3Vector PixelWorldNormal(const R2Point& image_position) const;
  R3Ray PixelWorldRay(const R2Point& image_position) const;
  RNScalar PixelChannelValue(const R2Point& image_position, int channel_index) const;

  // Camera property functions
  R4Matrix ProjectionMatrix(RNLength neardist = 0.1, RNLength fardist = 100.0) const;
  R4Matrix ModelViewMatrix(void) const;
  R3Point WorldViewpoint(void) const;
  R3Vector WorldTowards(void) const;
  R3Vector WorldRight(void) const;
  R3Vector WorldUp(void) const;
  R3Box WorldBBox(void) const;
  RNScalar Timestamp(void) const;
  RNAngle XFov(void) const;
  RNAngle YFov(void) const;
  
  // Transformation property functions
  const R3Affine& CameraToWorld(void) const;
  const R4Matrix Extrinsics(void) const;
  const R3Matrix& Intrinsics(void) const;

  // Manipulation functions
  virtual void SetNPixels(int nx, int ny);
  virtual void SetPixelColor(int ix, int iy, const RNRgb& color);
  virtual void SetPixelDepth(int ix, int iy, RNScalar depth);
  virtual void SetPixelCategory(int ix, int iy, RNScalar category);
  virtual void SetPixelInstance(int ix, int iy, RNScalar instance);
  virtual void SetPixelChannelValue(int ix, int iy, int channel_index, RNScalar value);
  virtual void SetChannel(int channel_index, const R2Grid& image);
  virtual void SetDepthChannel(const R2Grid& depth_image);
  virtual void SetColorChannels(const R2Image& color_image);
  virtual void SetCategoryChannel(const R2Grid& category_image);
  virtual void SetInstanceChannel(const R2Grid& instance_image);
  virtual void SetCameraToWorld(const R3Affine& transformation);
  virtual void SetExtrinsics(const R4Matrix& matrix);
  virtual void SetIntrinsics(const R3Matrix& matrix);
  virtual void SetTimestamp(RNScalar timestamp);
  virtual void Transform(const R3Transformation& transformation);

  // Transformation functions
  virtual int TransformImageToCamera(const R2Point& image_position, R3Point& camera_position) const;
  virtual int TransformCameraToImage(const R3Point& camera_position, R2Point& image_position) const;
  virtual int TransformCameraToWorld(const R3Point& camera_position, R3Point& world_position) const;
  virtual int TransformWorldToCamera(const R3Point& world_position, R3Point& camera_position) const;
  virtual int TransformWorldToImage(const R3Point& world_position, R2Point& image_position) const;

  // Draw functions
  virtual void Draw(int color_scheme = RGBD_PHOTO_COLOR_SCHEME) const;
  virtual void DrawCamera(int color_scheme = RGBD_PHOTO_COLOR_SCHEME, RNLength radius = 0.25) const;
  virtual void DrawBBox(int color_scheme = RGBD_INDEX_COLOR_SCHEME) const;
  virtual void DrawImage(int color_scheme = RGBD_PHOTO_COLOR_SCHEME, RNLength depth = 0) const;
  virtual void DrawPoints(int color_scheme = RGBD_PHOTO_COLOR_SCHEME, int skip = 1) const;
  virtual void DrawSurfels(int color_scheme = RGBD_PHOTO_COLOR_SCHEME, int skip = 1) const;
  virtual void DrawQuads(int color_scheme = RGBD_PHOTO_COLOR_SCHEME, int skip = 1) const;
  virtual void LoadColor(int color_scheme = RGBD_PHOTO_COLOR_SCHEME) const;

  // Create/read/write/release functions
  virtual int ReadChannels(void);
  virtual int ReadColorChannels(void);
  virtual int ReadDepthChannel(void);
  virtual int ReadCategoryChannel(void);
  virtual int ReadInstanceChannel(void);
  virtual int WriteChannels(void);
  virtual int WriteColorChannels(void);
  virtual int WriteDepthChannel(void);
  virtual int WriteCategoryChannel(void);
  virtual int WriteInstanceChannel(void);
  virtual int ReleaseChannels(void);
  virtual int ReleaseColorChannels(void);
  virtual int ReleaseDepthChannel(void);
  virtual int ReleaseCategoryChannel(void);
  virtual int ReleaseInstanceChannel(void);

public:
  // Channel access
  R2Image ColorChannels(void) const;

  // Channel creation
  virtual int CreateChannel(int channel, const R2Grid& image);
  virtual int CreateColorChannels(const R2Image& color_image);
  virtual int CreateDepthChannel(const R2Grid& image);
  virtual int CreateCategoryChannel(const R2Grid& image);
  virtual int CreateInstanceChannel(const R2Grid& image);
  
  // Filename access functions
  const char *Name(void) const;
  const char *ColorFilename(void) const;
  const char *DepthFilename(void) const;
  const char *CategoryFilename(void) const;
  const char *InstanceFilename(void) const;

  // Filename manipulation functions (do not read files)
  virtual void SetName(const char *name);
  virtual void SetColorFilename(const char *filename);
  virtual void SetDepthFilename(const char *filename);
  virtual void SetCategoryFilename(const char *filename);
  virtual void SetInstanceFilename(const char *filename);

  // Update functions
  virtual void InvalidateWorldBBox(void);
  virtual void InvalidateOpenGL(void);
  virtual void UpdateOpenGL(void);

  // Other functions
  int ComputeMesh(R3Mesh& mesh, RNLength max_depth = 0, RNLength max_silhouette_factor = 0.1);

  // User data functions
  void *UserData(void) const;
  void SetUserData(void *data);

private:
  // Internal variables
  friend class RGBDConfiguration;
  RGBDConfiguration *configuration;
  int configuration_index;
  RNArray<R2Grid *> channels;
  int width, height;
  R3Affine camera_to_world;
  R3Matrix intrinsics;
  R3Box world_bbox;
  RNScalar timestamp;
  int opengl_texture_id;
  char *name;
  char *color_filename;
  char *depth_filename;
  char *category_filename;
  char *instance_filename;
  int color_resident_count;
  int depth_resident_count;
  int category_resident_count;
  int instance_resident_count;
  void *data;
};



////////////////////////////////////////////////////////////////////////
// Inline functions
////////////////////////////////////////////////////////////////////////

inline RGBDConfiguration *RGBDImage::
Configuration(void) const
{
  // Return configuration of which this image is part
  return configuration;
}


inline int RGBDImage::
ConfigurationIndex(void) const
{
  // Return index of this image in configuration (i.e., where configuration->Image(k) == this)
  return configuration_index;
}



inline int RGBDImage::
NChannels(void) const
{
  // Return number of channels
  return channels.NEntries();
}



inline R2Grid *RGBDImage::
Channel(int channel_index) const
{
  // Check channel
  if (channel_index >= channels.NEntries()) return NULL;
  
  // Check channel
  if (!channels[channel_index]) return NULL;
  
  // Return channel
  return channels[channel_index];
}



inline R2Grid *RGBDImage::
RedChannel(void) const
{
  // Return channel
  return Channel(RGBD_RED_CHANNEL);
}



inline R2Grid *RGBDImage::
GreenChannel(void) const
{
  // Return channel
  return Channel(RGBD_GREEN_CHANNEL);
}



inline R2Grid *RGBDImage::
BlueChannel(void) const
{
  // Return channel
  return Channel(RGBD_BLUE_CHANNEL);
}



inline R2Grid *RGBDImage::
DepthChannel(void) const
{
  // Return channel
  return Channel(RGBD_DEPTH_CHANNEL);
}



inline R2Grid *RGBDImage::
CategoryChannel(void) const
{
  // Return channel
  return Channel(RGBD_CATEGORY_CHANNEL);
}



inline R2Grid *RGBDImage::
InstanceChannel(void) const
{
  // Return channel
  return Channel(RGBD_INSTANCE_CHANNEL);
}



inline int RGBDImage::
NPixels(int dim) const
{
  // Return number of pixels in dimension
  assert((dim >= 0) && (dim < 2));
  return (dim == 0) ? width : height;
}



inline RNScalar RGBDImage::
PixelChannelValue(int ix, int iy, int channel_index) const
{
  // Check channel
  if (!channels[channel_index]) {
    RNFail("RGBD Channel is not resident in memory -- cannot get value\n");
    return R2_GRID_UNKNOWN_VALUE;
  }
  
  // Return value in channel at pixel position
  return channels[channel_index]->GridValue(ix, iy);
}



inline RNScalar RGBDImage::
PixelDepth(int ix, int iy) const
{
  // Return depth of pixel
  if (NChannels() <= RGBD_DEPTH_CHANNEL) return 0;
  return PixelChannelValue(ix, iy, RGBD_DEPTH_CHANNEL);
}



inline RNScalar RGBDImage::
PixelCategory(int ix, int iy) const
{
  // Return category of pixel
  if (NChannels() <= RGBD_CATEGORY_CHANNEL) return 0;
  return PixelChannelValue(ix, iy, RGBD_CATEGORY_CHANNEL);
}



inline RNScalar RGBDImage::
PixelInstance(int ix, int iy) const
{
  // Return instance of pixel
  if (NChannels() <= RGBD_INSTANCE_CHANNEL) return 0;
  return PixelChannelValue(ix, iy, RGBD_INSTANCE_CHANNEL);
}



inline RNScalar RGBDImage::
PixelChannelValue(const R2Point& image_position, int channel_index) const
{
  // Return value in channel at pixel position
  return channels[channel_index]->GridValue(image_position);
}



inline RNScalar RGBDImage::
PixelDepth(const R2Point& image_position) const
{
  // Return depth of pixel
  if (NChannels() < 4) return 0;
  return PixelChannelValue(image_position, RGBD_DEPTH_CHANNEL);
}



inline R3Point RGBDImage::
PixelWorldPosition(int ix, int iy) const
{
  // Return position of pixel in world coordinates
  R2Point image_position(ix, iy);
  return PixelWorldPosition(image_position);
}



inline R3Vector RGBDImage::
PixelWorldNormal(const R2Point& image_position) const
{
  // Return position of pixel in world coordinates
  return PixelWorldNormal((int) image_position.X(), (int) image_position.Y());
}



inline R3Point RGBDImage::
WorldViewpoint(void) const
{
  // Return viewpoint of camera
  return camera_to_world.Matrix() * R3zero_point;
}



inline R3Vector RGBDImage::
WorldTowards(void) const
{
  // Return towards vector of camera
  return camera_to_world.Matrix() * R3negz_vector;
}



inline R3Vector RGBDImage::
WorldRight(void) const
{
  // Return right vector of camera
  return camera_to_world.Matrix() * R3posx_vector;
}



inline R3Vector RGBDImage::
WorldUp(void) const
{
  // Return up vector of camera
  return camera_to_world.Matrix() * R3posy_vector;
}



inline RNScalar RGBDImage::
Timestamp(void) const
{
  // Return timestamp
  return timestamp;
}



inline RNAngle RGBDImage::
XFov(void) const
{
  // Return field of view (half angle) in horizontal direction
  RNScalar xfocal = Intrinsics()[0][0];
  if (RNIsZero(xfocal)) return RN_INFINITY;
  else return atan(0.5 * NPixels(RN_X) / xfocal);
}



inline RNAngle RGBDImage::
YFov(void) const
{
  // Return field of view (half angle) in vertical direction
  RNScalar yfocal = Intrinsics()[1][1];
  if (RNIsZero(yfocal)) return RN_INFINITY;
  else return atan(0.5 * NPixels(RN_Y) / yfocal);
}



inline const R3Affine& RGBDImage::
CameraToWorld(void) const
{
  // Return matrix that transforms camera coordinates to world coordinates
  return camera_to_world;
}



inline const R4Matrix RGBDImage::
Extrinsics(void) const
{
  // Return extrinsics matrix (transforms camera coordinates to world coordinates)
  RNAbort("Not implemented for now");
  return camera_to_world.Matrix().Inverse();
}



inline const R3Matrix& RGBDImage::
Intrinsics(void) const
{
  // Return intrinsics matrix
  return intrinsics;
}



inline const char *RGBDImage::
Name(void) const
{
  // Return name
  return name;
}



inline const char *RGBDImage::
ColorFilename(void) const
{
  // Return color filename
  return color_filename;
}



inline const char *RGBDImage::
DepthFilename(void) const
{
  // Return depth filename
  return depth_filename;
}



inline const char *RGBDImage::
CategoryFilename(void) const
{
  // Return category filename
  return category_filename;
}



inline const char *RGBDImage::
InstanceFilename(void) const
{
  // Return instance filename
  return instance_filename;
}



inline void *RGBDImage::
UserData(void) const
{
  // Return user data
  return data;
}



inline void RGBDImage::
SetUserData(void *data)
{
  // Set user data
  this->data = data;
}



// End namespace
}


// End include guard
#endif
