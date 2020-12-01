/* Include file for the R3 surfel class */
#ifndef __R3__SURFEL__H__
#define __R3__SURFEL__H__



/* Begin namespace */
namespace gaps {



/* Class definition */

class R3Surfel {
public:
  // Constructor functions
  R3Surfel(void);
  R3Surfel(float px, float py, float pz, 
    unsigned char r = 0, unsigned char g = 0, unsigned char b = 0, 
    RNBoolean aerial = FALSE);
  R3Surfel(float px, float py, float pz, 
    float nx, float ny, float nz, 
    float radius = 0, 
    unsigned char r = 0, unsigned char g = 0, unsigned char b = 0, 
    unsigned char flags = 0);
  R3Surfel(float px, float py, float pz, 
    float nx, float ny, float nz, 
    float tx, float ty, float tz, 
    float radius0, float radius1,
    unsigned char r, unsigned char g, unsigned char b, 
    float timestamp, unsigned int identifier,
    unsigned char flags = 0);

  // Position property functions
  // NOTE THAT THESE COORDINATES ARE RELATIVE TO THE BLOCK ORIGIN
  // TO GET THE WORLD COORDINATES, YOU MUST ADD THE BLOCK ORIGIN
  float PX(void) const;
  float PY(void) const;
  float PZ(void) const;
  float PositionCoord(int dimension) const;

  // Normal property functions
  float NX(void) const;
  float NY(void) const;
  float NZ(void) const;
  float NormalCoord(int dimension) const;

  // Tangent property functions
  float TX(void) const;
  float TY(void) const;
  float TZ(void) const;
  float TangentCoord(int dimension) const;

  // Color property functions
  unsigned char R(void) const;
  unsigned char G(void) const;
  unsigned char B(void) const;
  RNRgb Rgb(void) const;

  // Radius property functions
  float Radius(void) const;
  float Radius(int axis) const;
  float AspectRatio(void) const;

  // Timestamp property functions
  float Timestamp(void) const;
  
  // Identifier property functions (can be used for anything)
  unsigned int Identifier(void) const;

  // Other property functions
  RNBoolean IsActive(void) const;
  RNBoolean IsMarked(void) const;
  RNBoolean IsAerial(void) const;
  RNBoolean IsTerrestrial(void) const;
  RNBoolean IsOriented(void) const;
  RNBoolean IsIsotropic(void) const;
  RNBoolean IsOnSilhouetteBoundary(void) const;
  RNBoolean IsOnShadowBoundary(void) const;
  RNBoolean IsOnBorderBoundary(void) const;
  RNBoolean IsOnBoundary(void) const;
  RNBoolean HasNormal(void) const;
  RNBoolean HasTangent(void) const;
  RNBoolean HasRadius(void) const;
  unsigned char Flags(void) const;

  // Manipulation functions
  void SetPosition(float x, float y, float z);
  void SetPosition(const float *xyz);
  void SetNormal(float x, float y, float z);
  void SetNormal(const float *xyz);
  void SetTangent(float x, float y, float z);
  void SetTangent(const float *xyz);
  void SetRadius(float radius);
  void SetRadius(int axis, float radius);
  void SetColor(unsigned char r, unsigned char g, unsigned char b);
  void SetColor(const unsigned char *rgb);
  void SetColor(const RNRgb& rgb);
  void SetTimestamp(float timestamp);
  void SetIdentifier(unsigned int identifier);
  void SetAerial(RNBoolean aerial = TRUE);
  void SetSilhouetteBoundary(RNBoolean boundary = TRUE);
  void SetShadowBoundary(RNBoolean boundary = TRUE);
  void SetBorderBoundary(RNBoolean boundary = TRUE);
  void SetMark(RNBoolean mark = TRUE);
  void SetActive(RNBoolean active = TRUE);
  void SetFlags(unsigned char flags);

  // Draw functions
  void Draw(RNFlags flags = R3_SURFEL_DEFAULT_DRAW_FLAGS) const;


////////////////////////////////////////////////////////////////////////
// INTERNAL STUFF BELOW HERE
////////////////////////////////////////////////////////////////////////

  // For backward compatibility
  float X(void) const;
  float Y(void) const;
  float Z(void) const;

  // For backward compatibility
  const float *Coords(void) const;
  float Coord(int dimension) const;

  // For backward compatibility
  void SetCoords(float x, float y, float z);
  void SetCoords(const float *xyz);

  // For backward compatibility
  const unsigned char *Color(void) const;

  // Do not use these 
  const float *PositionPtr(void) const;
  const RNInt16 *NormalPtr(void) const;
  const RNInt16 *TangentPtr(void) const;
  const RNUInt16 *RadiusPtr(void) const;
  const unsigned char *ColorPtr(void) const;

private:
  // Internal data
  RNScalar32 position[3];
  RNScalar32 timestamp;
  RNInt16 normal[3]; // x 2^15-1 (32767)
  RNInt16 tangent[3]; // x 2^15-1 (32767)
  RNUInt16 radius[2]; // x 2^13 (8192)
  RNUInt32 identifier;
  RNUChar8 color[3];
  RNUChar8 flags;

  // Friend database so that it can swap endian when reading
  friend class R3SurfelDatabase;
};



////////////////////////////////////////////////////////////////////////
// SURFEL FLAGS
////////////////////////////////////////////////////////////////////////

#define R3_SURFEL_ACTIVE_FLAG               0x04
#define R3_SURFEL_AERIAL_FLAG               0x02
#define R3_SURFEL_TANGENT_FLAG              0x01
#define R3_SURFEL_MARKED_FLAG               0x80
#define R3_SURFEL_NORMAL_FLAG               0x40
#define R3_SURFEL_SHADOW_BOUNDARY_FLAG      0x20
#define R3_SURFEL_SILHOUETTE_BOUNDARY_FLAG  0x10
#define R3_SURFEL_BORDER_BOUNDARY_FLAG      0x30
#define R3_SURFEL_BOUNDARY_FLAGS            0x30



////////////////////////////////////////////////////////////////////////
// INLINE FUNCTION DEFINITIONS
////////////////////////////////////////////////////////////////////////

inline float R3Surfel::
PX(void) const
{
  // Return X coordinate
  return position[0];
}



inline float R3Surfel::
PY(void) const
{
  // Return Y coordinate
  return position[1];
}



inline float R3Surfel::
PZ(void) const
{
  // Return Z coordinate
  return position[2];
}



inline float R3Surfel::
PositionCoord(int dimension) const
{
  // Return coordinate
  return position[dimension];
}



inline float R3Surfel::
NX(void) const
{
  // Return X normal coordinate
  return (1.0F/32767.0F) * normal[0];
}



inline float R3Surfel::
NY(void) const
{
  // Return Y normal coordinate
  return (1.0F/32767.0F) * normal[1];
}



inline float R3Surfel::
NZ(void) const
{
  // Return Z normal coordinate
  return (1.0F/32767.0F) * normal[2];
}



inline float R3Surfel::
NormalCoord(int dimension) const
{
  // Return normal coordinate
  return (1.0F/32767.0F) * normal[dimension];
}



inline float R3Surfel::
TX(void) const
{
  // Return X tangent coordinate
  return (1.0F/32767.0F) * tangent[0];
}



inline float R3Surfel::
TY(void) const
{
  // Return Y tangent coordinate
  return (1.0F/32767.0F) * tangent[1];
}



inline float R3Surfel::
TZ(void) const
{
  // Return Z tangent coordinate
  return (1.0F/32767.0F) * tangent[2];
}



inline float R3Surfel::
TangentCoord(int dimension) const
{
  // Return tangent coordinate
  return (1.0F/32767.0F) * tangent[dimension];
}



inline float R3Surfel::
Radius(void) const
{
  // Return longer radius
  return (1.0F/8192.0F) * radius[0];
}



inline float R3Surfel::
Radius(int axis) const
{
  // Return radius along axis
  assert((axis >= 0) && (axis <= 1));
  return (1.0F/8192.0F) * radius[axis];
}



inline float R3Surfel::
AspectRatio(void) const
{
  // Return radius[0] / radius[1]
  if (radius[1] == 0) return 1.0F;
  return (float) radius[0] / (float) radius[1];
}



inline unsigned char R3Surfel::
R(void) const
{
  // Return red component of color
  return color[0];
}



inline unsigned char R3Surfel::
G(void) const
{
  // Return green component of color
  return color[1];
}



inline unsigned char R3Surfel::
B(void) const
{
  // Return blue component of color
  return color[2];
}



inline const unsigned char *R3Surfel::
Color(void) const
{
  // Return pointer to color
  return color;
}



inline RNRgb R3Surfel::
Rgb(void) const
{
  // Return RGB
  const double scale = 1.0 / 255.0;
  return RNRgb(scale*color[0], scale*color[1], scale*color[2]);
}



inline float R3Surfel::
Timestamp(void) const
{
  // Return timestamp
  return timestamp;
}



inline unsigned int R3Surfel::
Identifier(void) const
{
  // Return identifier
  return identifier;
}



inline RNBoolean R3Surfel::
IsActive(void) const
{
  // Return whether point is active (not previously deleted or subsumed by some other point)
  return flags & R3_SURFEL_ACTIVE_FLAG;
}



inline RNBoolean R3Surfel::
IsAerial(void) const
{
  // Return whether point was captured with aerial scanner
  return flags & R3_SURFEL_AERIAL_FLAG;
}



inline RNBoolean R3Surfel::
IsTerrestrial(void) const
{
  // Return whether point was captured with terrestrial scanner
  return (!IsAerial());
}



inline RNBoolean R3Surfel::
IsOriented(void) const
{
  // Return whether point is oriented
  return (HasNormal() & HasTangent());
}



inline RNBoolean R3Surfel::
IsIsotropic(void) const
{
  // Return whether point is isotropic
  if ((radius[1] == 0) || (radius[0] == radius[1])) return TRUE;
  else return FALSE;
}



inline RNBoolean R3Surfel::
IsOnSilhouetteBoundary(void) const
{
  // Return whether point is on a silhouette boundary
  return ((flags & R3_SURFEL_BOUNDARY_FLAGS) == R3_SURFEL_SILHOUETTE_BOUNDARY_FLAG);
}



inline RNBoolean R3Surfel::
IsOnShadowBoundary(void) const
{
  // Return whether point is on a shadow boundary
  return ((flags & R3_SURFEL_BOUNDARY_FLAGS) == R3_SURFEL_SHADOW_BOUNDARY_FLAG);
}



inline RNBoolean R3Surfel::
IsOnBorderBoundary(void) const
{
  // Return whether point is on a shadow boundary
  return ((flags & R3_SURFEL_BOUNDARY_FLAGS) == R3_SURFEL_BORDER_BOUNDARY_FLAG);
}



inline RNBoolean R3Surfel::
IsOnBoundary(void) const
{
  // Return whether point is on a boundary
  return (flags & R3_SURFEL_BOUNDARY_FLAGS);
}



inline RNBoolean R3Surfel::
IsMarked(void) const
{
  // Return whether surfel is marked (useful for set and traversal operations)
  return flags & R3_SURFEL_MARKED_FLAG;
}



inline RNBoolean R3Surfel::
HasNormal(void) const
{
  // Return whether surfel has normal
  return flags & R3_SURFEL_NORMAL_FLAG;
}



inline RNBoolean R3Surfel::
HasTangent(void) const
{
  // Return whether surfel has tangent
  return flags & R3_SURFEL_TANGENT_FLAG;
}



inline RNBoolean R3Surfel::
HasRadius(void) const
{
  // Return whether surfel has non-zero radius
  return (radius[0] > 0) ? 1 : 0;
}



inline unsigned char R3Surfel::
Flags(void) const
{
  // Return bit-encoded status flags
  return flags;
}



inline void R3Surfel::
SetPosition(float x, float y, float z)
{
  // Set position
  position[0] = x;
  position[1] = y;
  position[2] = z;
}



inline void R3Surfel::
SetPosition(const float *xyz)
{
  // Set position
  position[0] = xyz[0];
  position[1] = xyz[1];
  position[2] = xyz[2];
}



inline void R3Surfel::
SetNormal(float x, float y, float z)
{
  // Set normal
  normal[0] = (RNInt16) (32767.0 * x + 0.5);
  normal[1] = (RNInt16) (32767.0 * y + 0.5);
  normal[2] = (RNInt16) (32767.0 * z + 0.5);

  // Update flags
  if ((x != 0) || (y != 0) || (z != 0)) flags |= R3_SURFEL_NORMAL_FLAG;
}



inline void R3Surfel::
SetNormal(const float *xyz)
{
  // Set normal
  SetNormal(xyz[0], xyz[1], xyz[2]);
}



inline void R3Surfel::
SetTangent(float x, float y, float z)
{
  // Set tangent
  tangent[0] = (RNInt16) (32767.0 * x + 0.5);
  tangent[1] = (RNInt16) (32767.0 * y + 0.5);
  tangent[2] = (RNInt16) (32767.0 * z + 0.5);

  // Update flags
  if ((x != 0) || (y != 0) || (z != 0)) flags |= R3_SURFEL_TANGENT_FLAG;
}



inline void R3Surfel::
SetTangent(const float *xyz)
{
  // Set tangent
  SetTangent(xyz[0], xyz[1], xyz[2]);
}



inline void R3Surfel::
SetRadius(float radius)
{
  // Set radius
  float r = 8192.0 * radius + 0.5;
  if (r > 65535) r = 65535;
  this->radius[0] = (RNUInt16) r;
  this->radius[1] = this->radius[0];
}



inline void R3Surfel::
SetRadius(int axis, float radius)
{
  // Set radius
  assert((axis >= 0) && (axis <= 1));
  float r = 8192.0 * radius + 0.5;
  if (r > 65535) r = 65535;
  this->radius[axis] = (RNUInt16) r;
}



inline void R3Surfel::
SetColor(unsigned char r, unsigned char g, unsigned char b)
{
  // Set color
  color[0] = r;
  color[1] = g;
  color[2] = b;
}



inline void R3Surfel::
SetColor(const unsigned char *rgb)
{
  // Set color
  color[0] = rgb[0];
  color[1] = rgb[1];
  color[2] = rgb[2];
}



inline void R3Surfel::
SetColor(const RNRgb& rgb)
{
  // Set color
  color[0] = (unsigned char) (255.0 * rgb.R());
  color[1] = (unsigned char) (255.0 * rgb.G());
  color[2] = (unsigned char) (255.0 * rgb.B());
}



inline void R3Surfel::
SetTimestamp(float timestamp)
{
  // Set timestamp
  this->timestamp = timestamp;
}



inline void R3Surfel::
SetIdentifier(unsigned int identifier)
{
  // Set identifier
  this->identifier = identifier;
}



inline void R3Surfel::
SetActive(RNBoolean active)
{
  // Set whether point is active
  if (active) flags |= R3_SURFEL_ACTIVE_FLAG;
  else flags &= ~R3_SURFEL_ACTIVE_FLAG;
}



inline void R3Surfel::
SetAerial(RNBoolean aerial)
{
  // Set whether point was captured with aerial scanner
  if (aerial) flags |= R3_SURFEL_AERIAL_FLAG;
  else flags &= ~R3_SURFEL_AERIAL_FLAG;
}



inline void R3Surfel::
SetSilhouetteBoundary(RNBoolean boundary)
{
  // Set whether point is on a silhouette boundary
  if (boundary) flags |= R3_SURFEL_SILHOUETTE_BOUNDARY_FLAG;
  else flags &= ~R3_SURFEL_SILHOUETTE_BOUNDARY_FLAG;
}



inline void R3Surfel::
SetShadowBoundary(RNBoolean boundary)
{
  // Set whether point is on a shadow boundary
  if (boundary) flags |= R3_SURFEL_SHADOW_BOUNDARY_FLAG;
  else flags &= ~R3_SURFEL_SHADOW_BOUNDARY_FLAG;
}



inline void R3Surfel::
SetBorderBoundary(RNBoolean boundary)
{
  // Set whether point is on a shadow boundary
  if (boundary) flags |= R3_SURFEL_BORDER_BOUNDARY_FLAG;
  else flags &= ~R3_SURFEL_BORDER_BOUNDARY_FLAG;
}



inline void R3Surfel::
SetMark(RNBoolean mark)
{
  // Set whether surfel is marked (useful for set and traversal operations)
  if (mark) flags |= R3_SURFEL_MARKED_FLAG;
  else flags &= ~R3_SURFEL_MARKED_FLAG;
}



inline void R3Surfel::
SetFlags(unsigned char flags)
{
  // Set flags -- only use this if you REALLY know what you are doing
  this->flags = flags;
}



inline void R3Surfel::
Draw(RNFlags flags) const
{
  // Draw point at surfel
  glBegin(GL_POINTS);
  if (flags[R3_SURFEL_COLOR_DRAW_FLAG]) glColor3ubv(color);
  glVertex3fv(position);
  glEnd();
}



////////////////////////////////////////////////////////////////////////
// Internal -- do not use
////////////////////////////////////////////////////////////////////////

inline const float *R3Surfel::
PositionPtr(void) const
{
  // Return pointer to position
  return position;
}



inline const RNInt16 *R3Surfel::
NormalPtr(void) const
{
  // Return pointer to normal
  return normal;
}



inline const RNInt16 *R3Surfel::
TangentPtr(void) const
{
  // Return pointer to tangent
  return tangent;
}



inline const RNUInt16 *R3Surfel::
RadiusPtr(void) const
{
  // Return pointer to radius
  return radius;
}



inline const unsigned char *R3Surfel::
ColorPtr(void) const
{
  // Return pointer to color
  return color;
}



////////////////////////////////////////////////////////////////////////
// For backwards compatibility -- do not use
////////////////////////////////////////////////////////////////////////

inline float R3Surfel::
X(void) const
{
  // Return X coordinate
  return PX();
}



inline float R3Surfel::
Y(void) const
{
  // Return Y coordinate
  return PY();
}



inline float R3Surfel::
Z(void) const
{
  // Return Z coordinate
  return PZ();
}



inline float R3Surfel::
Coord(int dimension) const
{
  // Return coordinate
  return PositionCoord(dimension);
}



inline const float *R3Surfel::
Coords(void) const
{
  // Return pointer to position
  return position;
}



inline void R3Surfel::
SetCoords(float x, float y, float z)
{
  // Set position
  SetPosition(x, y, z);
}



inline void R3Surfel::
SetCoords(const float *xyz)
{
  // Set position
  SetPosition(xyz);
}



// End namespace
}


// End include guard
#endif
