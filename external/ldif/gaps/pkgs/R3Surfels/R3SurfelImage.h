/* Include file for the R3 surfel image class */
#ifndef __R3__SURFEL__IMAGE__H__
#define __R3__SURFEL__IMAGE__H__



////////////////////////////////////////////////////////////////////////
// NAMESPACE 
////////////////////////////////////////////////////////////////////////

namespace gaps {



////////////////////////////////////////////////////////////////////////
// CLASS DEFINITION
////////////////////////////////////////////////////////////////////////

class R3SurfelImage {
public:
  //////////////////////////////////////////
  //// CONSTRUCTOR/DESTRUCTOR FUNCTIONS ////
  //////////////////////////////////////////

  // Constructor functions
  R3SurfelImage(const char *name = NULL);
  R3SurfelImage(const R3SurfelImage& image);

  // Destructor function
  virtual ~R3SurfelImage(void);


  ////////////////////////////
  //// PROPERTY FUNCTIONS ////
  ////////////////////////////

  // Pixel property access functions
  RNRgb PixelColor(int ix, int iy) const;
  RNScalar PixelRed(int ix, int iy) const;
  RNScalar PixelGreen(int ix, int iy) const;
  RNScalar PixelBlue(int ix, int iy) const;
  RNScalar PixelDepth(int ix, int iy) const;
  RNScalar PixelChannelValue(int ix, int iy, int channel_index) const;

  // Point property access functions
  RNRgb PixelColor(const R2Point& image_position) const;
  RNScalar PixelRed(const R2Point& image_position) const;
  RNScalar PixelGreen(const R2Point& image_position) const;
  RNScalar PixelBlue(const R2Point& image_position) const;
  R3Ray PixelWorldRay(const R2Point& image_position) const;
  RNScalar PixelChannelValue(const R2Point& image_position, int channel_index) const;

  // Channel access functions
  int NChannels(void) const;
  const R2Grid *Channel(int channel_index) const;
  const R2Grid *RedChannel(void) const;
  const R2Grid *GreenChannel(void) const;
  const R2Grid *BlueChannel(void) const;
  const R2Grid *DepthChannel(void) const;
  R2Image ColorChannels(void) const;
  
  // Camera intrinsics functions
  int ImageWidth(void) const;
  int ImageHeight(void) const;
  const R2Point& ImageCenter(void) const;
  RNLength XFocal(void) const;
  RNLength YFocal(void) const;
  RNAngle XFOV(void) const;
  RNAngle YFOV(void) const;
  R3Matrix Intrinsics(void) const;
  R4Matrix ProjectionMatrix(RNScalar neardist, RNScalar fardist) const;
  
  // Camera extrinsics functions
  const R3CoordSystem& Pose(void) const;
  const R3Point& Viewpoint(void) const;
  R3Vector Towards(void) const;
  const R3Vector& Up(void) const;
  const R3Vector& Right(void) const;
  R4Matrix CameraToWorld(void) const;
  R4Matrix Extrinsics(void) const;

  // Timestamp property functions
  RNScalar Timestamp(void) const;

  // Name property functions
  const char *Name(void) const;

  // User data property functions
  RNFlags Flags(void) const;
  void *Data(void) const;


  //////////////////////////
  //// ACCESS FUNCTIONS ////
  //////////////////////////

  // Scene access functions
  R3SurfelScene *Scene(void) const;
  int SceneIndex(void) const;

  // Scan access functions (can be NULL)
  R3SurfelScan *Scan(void) const;


  /////////////////////////////////////////
  //// PROPERTY MANIPULATION FUNCTIONS ////
  /////////////////////////////////////////

  // Channel manipulation functions
  virtual void SetChannel(int channel_index, const R2Grid& channel);
  virtual void SetRedChannel(const R2Grid& channel);
  virtual void SetGreenChannel(const R2Grid& channel);
  virtual void SetBlueChannel(const R2Grid& channel);
  virtual void SetDepthChannel(const R2Grid& channel);
  virtual void SetColorChannels(const R2Image& image);
  virtual void RemoveChannel(int channel_index);
  
  // Pose manipulation functions
  virtual void SetPose(const R3CoordSystem& pose);
  virtual void SetViewpoint(const R3Point& viewpoint);
  virtual void SetOrientation(const R3Vector& towards, const R3Vector& up);
  virtual void SetFocalLengths(RNLength focal_length);
  virtual void SetXFocal(RNLength focal_length);
  virtual void SetYFocal(RNLength focal_length);

  // Timestamp manipulation functions
  virtual void SetTimestamp(RNScalar timestamp);

  // Name manipulation functions
  virtual void SetName(const char *name);

  // Image metadata manipulation functions
  virtual void SetImageDimensions(int width, int height);
  virtual void SetImageCenter(const R2Point& center);

  // Set scan (if surfels captured with image)
  virtual void SetScan(R3SurfelScan *scan);
  
  // User data manipulation functions
  virtual void SetFlags(RNFlags flags);
  virtual void SetData(void *data);


  /////////////////////////////////////////////
  //// COORDINATE TRANSFORMATION FUNCTIONS ////
  /////////////////////////////////////////////

  // Transform between coordinate systems
  R3Point TransformFromWorldToCamera(const R3Point& world_position) const;
  R2Point TransformFromWorldToImage(const R3Point& world_position) const;
  R3Point TransformFromCameraToWorld(const R3Point& camera_position) const;
  R2Point TransformFromCameraToImage(const R3Point& camera_position) const;
  R3Point TransformFromImageToWorld(const R2Point& image_position) const;
  R3Point TransformFromImageToCamera(const R2Point& image_position, RNLength depth = -1) const;

  // Check if pixel coordinates are within image bounds
  RNBoolean ContainsImagePosition(const R2Point& image_position) const;

  
  ///////////////////////////
  //// DISPLAY FUNCTIONS ////
  ///////////////////////////

  // Draw function
  virtual void Draw(RNFlags flags = R3_SURFEL_DEFAULT_DRAW_FLAGS) const;

  // Print function
  virtual void Print(FILE *fp = NULL, const char *prefix = NULL, const char *suffix = NULL) const;

  // Render image by projecting surfels into image
  int RenderImage(R2Image *color_image = NULL,
    R2Grid *depth_image = NULL, R2Grid *height_image = NULL,
    R2Grid *xnormal_image = NULL, R2Grid *ynormal_image = NULL, R2Grid *znormal_image = NULL,
    R2Grid *label_image = NULL, R2Grid *object_image = NULL,
    R2Grid *node_image = NULL, R2Grid *block_image = NULL) const;

  
  ////////////////////////////////////////////////////////////////////////
  // INTERNAL STUFF BELOW HERE
  ////////////////////////////////////////////////////////////////////////

  // For backward compatibility
  R2Point ImagePosition(const R3Point& world_position) const;

protected:
  // Internal data
  friend class R3SurfelScene;
  friend class R3SurfelScan;
  R3SurfelScene *scene;
  int scene_index;
  R3SurfelScan *scan;
  RNArray<R2Grid *> channels;
  R3CoordSystem pose;
  RNScalar timestamp;
  int image_width, image_height;
  R2Point image_center;
  RNLength xfocal, yfocal;
  char *name;
  RNFlags flags;
  void *data;
};



////////////////////////////////////////////////////////////////////////
// Channel indices
////////////////////////////////////////////////////////////////////////

enum {
  R3_SURFEL_RED_CHANNEL,
  R3_SURFEL_GREEN_CHANNEL,
  R3_SURFEL_BLUE_CHANNEL,
  R3_SURFEL_DEPTH_CHANNEL,
  R3_SURFEL_USER_CHANNEL,
  R3_SURFEL_NUM_CHANNELS
};


  
////////////////////////////////////////////////////////////////////////
// INLINE FUNCTION DEFINITIONS
////////////////////////////////////////////////////////////////////////

inline int R3SurfelImage::
NChannels(void) const
{
  // Return the number of channels (including NULL pointers)
  return channels.NEntries();
}

  

inline const R2Grid *R3SurfelImage::
Channel(int channel_index) const
{
  // Return the channel (may be NULL)
  if (channel_index < 0) return NULL;
  if (channel_index >= channels.NEntries()) return NULL;
  return channels[channel_index];
}

  

inline const R2Grid *R3SurfelImage::
RedChannel(void) const
{
  // Return the red channel (may be NULL)
  return Channel(R3_SURFEL_RED_CHANNEL);
}



inline const R2Grid *R3SurfelImage::
GreenChannel(void) const
{
  // Return the green channel (may be NULL)
  return Channel(R3_SURFEL_GREEN_CHANNEL);
}



inline const R2Grid *R3SurfelImage::
BlueChannel(void) const
{
  // Return the blue channel (may be NULL)
  return Channel(R3_SURFEL_BLUE_CHANNEL);
}



inline const R2Grid *R3SurfelImage::
DepthChannel(void) const
{
  // Return the depth channel (may be NULL)
  return Channel(R3_SURFEL_DEPTH_CHANNEL);
}



inline RNScalar R3SurfelImage::
PixelRed(int ix, int iy) const
{
  // Return the red component of the pixel
  const R2Grid *red_channel = RedChannel();
  if (!red_channel) return -1;
  return red_channel->GridValue(ix, iy);
}



inline RNScalar R3SurfelImage::
PixelGreen(int ix, int iy) const
{
  // Return the green component of the pixel
  const R2Grid *green_channel = GreenChannel();
  if (!green_channel) return -1;
  return green_channel->GridValue(ix, iy);
}



inline RNScalar R3SurfelImage::
PixelBlue(int ix, int iy) const
{
  // Return the blue component of the pixel
  const R2Grid *blue_channel = BlueChannel();
  if (!blue_channel) return -1;
  return blue_channel->GridValue(ix, iy);
}



inline RNRgb R3SurfelImage::
PixelColor(int ix, int iy) const
{
  // Return the color of the pixel
  return RNRgb(PixelRed(ix,iy), PixelGreen(ix,iy), PixelBlue(ix,iy));
}



inline RNScalar R3SurfelImage::
PixelDepth(int ix, int iy) const
{
  // Return the depth of the pixel
  const R2Grid *depth_channel = DepthChannel();
  if (!depth_channel) return -1;
  return depth_channel->GridValue(ix, iy);
}



inline RNScalar R3SurfelImage::
PixelChannelValue(int ix, int iy, int channel_index) const
{
  // Return the channel value of the pixel
  const R2Grid *channel = Channel(channel_index);
  if (!channel) return -1;
  return channel->GridValue(ix, iy);
}



inline RNScalar R3SurfelImage::
PixelRed(const R2Point& image_position) const
{
  // Return the red component of the pixel
  const R2Grid *red_channel = RedChannel();
  if (!red_channel) return -1;
  return red_channel->GridValue(image_position);
}



inline RNScalar R3SurfelImage::
PixelGreen(const R2Point& image_position) const
{
  // Return the green component of the pixel
  const R2Grid *green_channel = GreenChannel();
  if (!green_channel) return -1;
  return green_channel->GridValue(image_position);
}



inline RNScalar R3SurfelImage::
PixelBlue(const R2Point& image_position) const
{
  // Return the blue component of the pixel
  const R2Grid *blue_channel = BlueChannel();
  if (!blue_channel) return -1;
  return blue_channel->GridValue(image_position);
}



inline RNRgb R3SurfelImage::
PixelColor(const R2Point& image_position) const
{
  // Return the color of the pixel
  return RNRgb(PixelRed(image_position), PixelGreen(image_position), PixelBlue(image_position));
}



inline RNScalar R3SurfelImage::
PixelChannelValue(const R2Point& image_position, int channel_index) const
{
  // Return the channel value of the pixel
  const R2Grid *channel = Channel(channel_index);
  if (!channel) return -1;
  return channel->GridValue(image_position);
}



inline RNLength R3SurfelImage::
XFocal(void) const
{
  // Return horizontal focal length in pixels
  return xfocal;
}



inline RNLength R3SurfelImage::
YFocal(void) const
{
  // Return vertical focal length in pixels
  return yfocal;
}



inline RNLength R3SurfelImage::
XFOV(void) const
{
  // Return half-angle for horizontal field of view
  if (xfocal <= 0) return 0.0;
  return atan(0.5*image_width/xfocal);
}



inline RNLength R3SurfelImage::
YFOV(void) const
{
  // Return half-angle for vertical field of view
  if (yfocal <= 0) return 0.0;
  return atan(0.5*image_height/yfocal);
}



inline int R3SurfelImage::
ImageWidth(void) const
{
  // Return image width
  return image_width;
}



inline int R3SurfelImage::
ImageHeight(void) const
{
  // Return image height
  return image_height;
}



inline const R2Point& R3SurfelImage::
ImageCenter(void) const
{
  // Return image center
  return image_center;
}



inline R3Matrix R3SurfelImage::
Intrinsics(void) const
{
  // Return intrinsics matrix
  return R3Matrix(xfocal, 0, image_center.X(),
                  0, yfocal, image_center.Y(),
                  0, 0, 1);
}



inline const R3CoordSystem& R3SurfelImage::
Pose(void) const
{
  // Return pose 
  return pose;
}



inline R4Matrix R3SurfelImage::
CameraToWorld(void) const
{
  // Return matrix going from camera coordinates to world coordinates
  return pose.Matrix();
}



inline R4Matrix R3SurfelImage::
Extrinsics(void) const
{
  // Return matrix going from world coordinates to camera coordinates
  return pose.InverseMatrix();
}



inline const R3Point& R3SurfelImage::
Viewpoint(void) const
{
  // Return pose viewpoint
  return pose.Origin();
}



inline R3Vector R3SurfelImage::
Towards(void) const
{
  // Return pose towards vector
  return -(pose.Axes().Axis(RN_Z));
}



inline const R3Vector& R3SurfelImage::
Up(void) const
{
  // Return pose up vector
  return pose.Axes().Axis(RN_Y);
}



inline const R3Vector& R3SurfelImage::
Right(void) const
{
  // Return pose right vector
  return pose.Axes().Axis(RN_X);
}



inline RNScalar R3SurfelImage::
Timestamp(void) const
{
  // Return timestamp
  return timestamp;
}



inline const char *R3SurfelImage::
Name(void) const
{
  // Return name
  return name;
}



inline RNFlags R3SurfelImage::
Flags(void) const
{
  // Return flags
  return flags;
}



inline void *R3SurfelImage::
Data(void) const
{
  // Return user data
  return data;
}



inline R3SurfelScene *R3SurfelImage::
Scene(void) const
{
  // Return scene this image is in
  return scene;
}



inline int R3SurfelImage::
SceneIndex(void) const
{
  // Return index in list of images associated with scene
  return scene_index;
}



inline R3SurfelScan *R3SurfelImage::
Scan(void) const
{
  // Return scan this image is associated (NULL if none)
  return scan;
}



inline R2Point R3SurfelImage::
ImagePosition(const R3Point& world_position) const
{
  // DO NOT USE -- only here for backward compatibility
  // Transform 3D point from world into image coordinate system
  return TransformFromWorldToImage(world_position);
}



// End namespace
}


// End include guard
#endif
