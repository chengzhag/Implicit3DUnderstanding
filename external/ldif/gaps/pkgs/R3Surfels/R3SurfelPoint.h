/* Include file for the R3 surfel point class */
#ifndef __R3__SURFEL__POINT__H__
#define __R3__SURFEL__POINT__H__



// Namespace 

namespace gaps {



/* Class definition */

class R3SurfelPoint {
public:
  // Constructor functions
  R3SurfelPoint(void);
  R3SurfelPoint(const R3SurfelPoint& surfel);
  R3SurfelPoint(R3SurfelBlock *block, const R3Surfel *surfel);
  R3SurfelPoint(R3SurfelBlock *block, int surfel_index);
  ~R3SurfelPoint(void);

  // Position property functions
  RNCoord PX(void) const;
  RNCoord PY(void) const;
  RNCoord PZ(void) const;
  RNCoord PositionCoord(int dimension) const;
  R3Point Position(void) const;

  // Normal property functions
  RNCoord NX(void) const;
  RNCoord NY(void) const;
  RNCoord NZ(void) const;
  RNCoord NormalCoord(int dimension) const;
  R3Vector Normal(void) const;

  // Tangent property functions
  RNCoord TX(void) const;
  RNCoord TY(void) const;
  RNCoord TZ(void) const;
  RNCoord TangentCoord(int dimension) const;
  R3Vector Tangent(void) const;

  // Color property functions
  unsigned char R(void) const;
  unsigned char G(void) const;
  unsigned char B(void) const;
  const unsigned char *Color(void) const;
  RNRgb Rgb(void) const;

  // Radius property functions
  float Radius(void) const;
  float Radius(int axis) const;

  // Timestamp property functions
  RNScalar Timestamp(void) const;

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

  // Access functions
  R3SurfelBlock *Block(void) const;
  int BlockIndex(void) const;
  const R3Surfel *Surfel(void) const;

  // Copy and reset functions
  void Copy(const R3SurfelPoint *point);
  void Reset(R3SurfelBlock *block, const R3Surfel *surfel);
  R3SurfelPoint& operator=(const R3SurfelPoint& point);

  // Manipulation functions
  void SetPosition(const R3Point& position);
  void SetNormal(const R3Vector& normal);
  void SetTangent(const R3Vector& tangent);
  void SetRadius(float radius);
  void SetRadius(int axis, float radius);
  void SetColor(const RNRgb& color);
  void SetTimestamp(RNScalar timestamp);
  void SetIdentifier(unsigned int identifier);
  void SetActive(RNBoolean active = TRUE);
  void SetAerial(RNBoolean aerial = TRUE);
  void SetMark(RNBoolean mark);

  // Draw functions
  void Draw(RNFlags flags = R3_SURFEL_DEFAULT_DRAW_FLAGS) const;

public:
  // Do not use these - for backwards compatibility
  // Position property functions
  RNCoord X(void) const;
  RNCoord Y(void) const;
  RNCoord Z(void) const;
  RNCoord Coord(int dimension) const;

private:
  R3SurfelBlock *block;
  const R3Surfel *surfel;
};



/* Public function declarations */

extern R3Point SurfelPointPosition(R3SurfelPoint *point, void *);



/* Inline functions */

inline RNCoord R3SurfelPoint::
PX(void) const
{
  // Return X coordinate of surfel point in global coordinate system
  assert(block && surfel);
 return surfel->PX() + block->PositionOrigin().X();
}



inline RNCoord R3SurfelPoint::
PY(void) const
{
  // Return Y coordinate of surfel point in global coordinate system
  assert(block && surfel);
  return surfel->PY() + block->PositionOrigin().Y();
}



inline RNCoord R3SurfelPoint::
PZ(void) const
{
  // Return Z coordinate of surfel point in global coordinate system
  assert(block && surfel);
  return surfel->PZ() + block->PositionOrigin().Z();
}



inline RNCoord R3SurfelPoint::
PositionCoord(int dimension) const
{
  // Return coordinate of surfel point in global coordinate system
  assert(block && surfel);
  return surfel->PositionCoord(dimension) + block->PositionOrigin().Coord(dimension);
}



inline R3Point R3SurfelPoint::
Position(void) const
{
  // Return position of surfel point in global coordinate system
  return R3Point(PX(), PY(), PZ());
}



inline RNCoord R3SurfelPoint::
NX(void) const
{
  // Return X normal of surfel 
  assert(block && surfel);
 return surfel->NX();
}



inline RNCoord R3SurfelPoint::
NY(void) const
{
  // Return Y normal of surfel 
  assert(block && surfel);
  return surfel->NY();
}



inline RNCoord R3SurfelPoint::
NZ(void) const
{
  // Return Z normal of surfel 
  assert(block && surfel);
  return surfel->NZ();
}



inline RNCoord R3SurfelPoint::
NormalCoord(int dimension) const
{
  // Return normal coordinate of surfel point 
  assert(block && surfel);
  return surfel->NormalCoord(dimension);
}



inline R3Vector R3SurfelPoint::
Normal(void) const
{
  // Return position of surfel point in global coordinate system
  return R3Vector(NX(), NY(), NZ());
}



inline RNCoord R3SurfelPoint::
TX(void) const
{
  // Return X tangent of surfel 
  assert(block && surfel);
 return surfel->TX();
}



inline RNCoord R3SurfelPoint::
TY(void) const
{
  // Return Y tangent of surfel 
  assert(block && surfel);
  return surfel->TY();
}



inline RNCoord R3SurfelPoint::
TZ(void) const
{
  // Return Z tangent of surfel 
  assert(block && surfel);
  return surfel->TZ();
}



inline RNCoord R3SurfelPoint::
TangentCoord(int dimension) const
{
  // Return tangent coordinate of surfel point 
  assert(block && surfel);
  return surfel->TangentCoord(dimension);
}



inline R3Vector R3SurfelPoint::
Tangent(void) const
{
  // Return position of surfel point in global coordinate system
  return R3Vector(TX(), TY(), TZ());
}



inline float R3SurfelPoint::
Radius(void) const
{
  // Return radius of surfel
  return surfel->Radius();
}



inline float R3SurfelPoint::
Radius(int axis) const
{
  // Return radius of surfel
  return surfel->Radius(axis);
}



inline unsigned char R3SurfelPoint::
R(void) const
{
  // Return red component of color
  return surfel->R();
}



inline unsigned char R3SurfelPoint::
G(void) const
{
  // Return green component of color
  return surfel->G();
}



inline unsigned char R3SurfelPoint::
B(void) const
{
  // Return blue component of color
  return surfel->B();
}



inline const unsigned char *R3SurfelPoint::
Color(void) const
{
  // Return pointer to color
  return surfel->Color();
}



inline RNRgb R3SurfelPoint::
Rgb(void) const
{
  // Return RGB
  return surfel->Rgb();
}



inline RNScalar R3SurfelPoint::
Timestamp(void) const
{
  // Return timestamp of surfel
  return surfel->Timestamp() + block->TimestampOrigin();
}



inline unsigned int R3SurfelPoint::
Identifier(void) const
{
  // Return identifier of surfel
  return surfel->Identifier();
}



inline RNBoolean R3SurfelPoint::
IsActive(void) const
{
  // Return whether point is active
  return surfel->IsActive();
}



inline RNBoolean R3SurfelPoint::
IsMarked(void) const
{
  // Return whether point is marked
  return surfel->IsMarked();
}



inline RNBoolean R3SurfelPoint::
IsAerial(void) const
{
  // Return whether point was captured with aerial scanner
  return surfel->IsAerial();
}



inline RNBoolean R3SurfelPoint::
IsTerrestrial(void) const
{
  // Return whether point was captured with terrestrial scanner
  return surfel->IsTerrestrial();
}



inline RNBoolean R3SurfelPoint::
IsOriented(void) const
{
  // Return whether point is oriented (has a normal and tangent)
  return surfel->IsOriented();
}



inline RNBoolean R3SurfelPoint::
IsIsotropic(void) const
{
  // Return whether point is isotropic (radius0 == radius1)
  return surfel->IsTerrestrial();
}



inline RNBoolean R3SurfelPoint::
IsOnSilhouetteBoundary(void) const
{
  // Return whether point is on silhouette boundary
  return surfel->IsOnSilhouetteBoundary();
}



inline RNBoolean R3SurfelPoint::
IsOnShadowBoundary(void) const
{
  // Return whether point is on shadow boundary
  return surfel->IsOnShadowBoundary();
}



inline RNBoolean R3SurfelPoint::
IsOnBorderBoundary(void) const
{
  // Return whether point is on border boundary
  return surfel->IsOnBorderBoundary();
}



inline RNBoolean R3SurfelPoint::
IsOnBoundary(void) const
{
  // Return whether point is on boundary
  return surfel->IsOnBoundary();
}



inline RNBoolean R3SurfelPoint::
HasNormal(void) const
{
  // Return whether point has normal
  return surfel->HasNormal();
}



inline RNBoolean R3SurfelPoint::
HasTangent(void) const
{
  // Return whether point has tangent
  return surfel->HasTangent();
}



inline RNBoolean R3SurfelPoint::
HasRadius(void) const
{
  // Return whether point has radius
  return surfel->HasRadius();
}



inline unsigned char R3SurfelPoint::
Flags(void) const
{
  // Return bit-encoded status flags
  return surfel->Flags();
}



inline R3SurfelBlock *R3SurfelPoint::
Block(void) const
{
  // Return block
  return block;
}



inline int R3SurfelPoint::
BlockIndex(void) const
{
  // Return index of surfel within block
  return surfel - block->Surfels();
}



inline const R3Surfel *R3SurfelPoint::
Surfel(void) const
{
  // Return surfel 
  return surfel;
}



////////////////////////////////////////////////////////////////////////
// Do not use these -- for backward compatibility
////////////////////////////////////////////////////////////////////////

inline RNCoord R3SurfelPoint::
X(void) const
{
  // Return X coordinate of surfel point in global coordinate system
  return PX();
}



inline RNCoord R3SurfelPoint::
Y(void) const
{
  // Return Y coordinate of surfel point in global coordinate system
  return PY();
}



inline RNCoord R3SurfelPoint::
Z(void) const
{
  // Return Z coordinate of surfel point in global coordinate system
  return PZ();
}



inline RNCoord R3SurfelPoint::
Coord(int dimension) const
{
  // Return coordinate of surfel point in global coordinate system
  return PositionCoord(dimension);
}



// End namespace
}


// End include guard
#endif
