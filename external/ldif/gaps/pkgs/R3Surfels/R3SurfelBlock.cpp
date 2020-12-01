/* Source file for the R3 surfel block class */



////////////////////////////////////////////////////////////////////////
// INCLUDE FILES
////////////////////////////////////////////////////////////////////////

#include "R3Surfels.h"



// Namespace

namespace gaps {



////////////////////////////////////////////////////////////////////////
// CONSTRUCTORS/DESTRUCTORS
////////////////////////////////////////////////////////////////////////

R3SurfelBlock::
R3SurfelBlock(void)
  : surfels(NULL),
    nsurfels(0),
    position_origin(0,0,0),
    bbox(FLT_MAX,FLT_MAX,FLT_MAX,-FLT_MAX,-FLT_MAX,-FLT_MAX),
    timestamp_origin(0),
    timestamp_range(FLT_MAX,-FLT_MAX),
    min_identifier(UINT_MAX),
    max_identifier(0),
    resolution(0),
    flags(0),
    data(NULL),
    database(NULL),
    database_index(-1),
    file_surfels_offset(0),
    file_surfels_count(0),
    file_read_count(0),
    node(NULL),
    opengl_id(0)
{
}



R3SurfelBlock::
R3SurfelBlock(int nsurfels)
  : surfels(NULL),
    nsurfels(nsurfels),
    position_origin(0,0,0),
    bbox(FLT_MAX,FLT_MAX,FLT_MAX,-FLT_MAX,-FLT_MAX,-FLT_MAX),
    timestamp_origin(0),
    timestamp_range(FLT_MAX,-FLT_MAX),
    min_identifier(UINT_MAX),
    max_identifier(0),
    resolution(0),
    flags(0),
    data(NULL),
    database(NULL),
    database_index(-1),
    file_surfels_offset(0),
    file_surfels_count(0),
    file_read_count(0),
    node(NULL),
    opengl_id(0)
{
  // Allocate surfels
  if (nsurfels > 0) {
    this->surfels = new R3Surfel [ nsurfels ];
  }
}



R3SurfelBlock::
R3SurfelBlock(const R3SurfelBlock& block)
  : surfels(NULL),
    nsurfels(block.nsurfels),
    position_origin(block.position_origin),
    bbox(block.bbox),
    timestamp_origin(block.timestamp_origin),
    timestamp_range(block.timestamp_range),
    min_identifier(block.min_identifier),
    max_identifier(block.max_identifier),
    resolution(block.resolution),
    flags(block.flags & R3_SURFEL_BLOCK_PROPERTY_FLAGS),
    data(NULL),
    database(NULL),
    database_index(-1),
    file_surfels_offset(0),
    file_surfels_count(0),
    file_read_count(0),
    node(NULL),
    opengl_id(0)
{
  // Copy surfels
  if (nsurfels > 0) {
    surfels = new R3Surfel [ nsurfels ];
    for (int i = 0; i < nsurfels; i++) {
      surfels[i] = block.surfels[i];
    }
  }
}



R3SurfelBlock::
R3SurfelBlock(const R3SurfelPointSet *set)
  : surfels(NULL),
    nsurfels(set->NPoints()),
    position_origin(set->Centroid()),
    bbox(set->BBox()),
    timestamp_origin(0),
    timestamp_range(FLT_MAX,-FLT_MAX),
    min_identifier(UINT_MAX),
    max_identifier(0),
    resolution(0),
    flags(R3_SURFEL_BLOCK_BBOX_UPTODATE_FLAG),
    data(NULL),
    database(NULL),
    database_index(-1),
    file_surfels_offset(0),
    file_surfels_count(0),
    file_read_count(0),
    node(NULL),
    opengl_id(0)
{
  // Copy points
  if (nsurfels > 0) {
    // Set timestamp range and origin
    this->timestamp_range = set->TimestampRange();
    if (!this->timestamp_range.IsEmpty()) {
      this->timestamp_origin = this->timestamp_range.Mid();
    }

    // Create surfels
    surfels = new R3Surfel [ nsurfels ];
    for (int i = 0; i < set->NPoints(); i++) {
      const R3SurfelPoint *point = set->Point(i);
      R3Point position = point->Position() - position_origin.Vector();
      R3Vector normal = point->Normal();
      R3Vector tangent = point->Tangent();
      R3Surfel *surfel = &surfels[i];
      surfel->SetPosition(position[0], position[1], position[2]);
      surfel->SetNormal(normal[0], normal[1], normal[2]);
      surfel->SetTangent(tangent[0], tangent[1], tangent[2]);
      surfel->SetColor(point->Color());
      surfel->SetRadius(0, point->Radius(0));
      surfel->SetRadius(1, point->Radius(1));
      surfel->SetTimestamp(point->Timestamp() - timestamp_origin);
      surfel->SetIdentifier(point->Identifier());
      surfel->SetFlags(point->Flags() & ~R3_SURFEL_MARKED_FLAG);
    }
  }
}



R3SurfelBlock::
R3SurfelBlock(const R3SurfelPointSet *set,
  const R3Point& position_origin, RNScalar timestamp_origin)
  : surfels(NULL),
    nsurfels(set->NPoints()),
    position_origin(position_origin),
    bbox(set->BBox()),
    timestamp_origin(0),
    timestamp_range(FLT_MAX,-FLT_MAX),
    min_identifier(UINT_MAX),
    max_identifier(0),
    resolution(0),
    flags(R3_SURFEL_BLOCK_BBOX_UPTODATE_FLAG),
    data(NULL),
    database(NULL),
    database_index(-1),
    file_surfels_offset(0),
    file_surfels_count(0),
    file_read_count(0),
    node(NULL),
    opengl_id(0)
{
  // Copy surfels
  if (nsurfels > 0) {
    // Set timestamp range and origin
    this->timestamp_range = set->TimestampRange();
    if (!this->timestamp_range.IsEmpty()) {
      this->timestamp_origin = this->timestamp_range.Mid();
    }

    // Copy surfels
    surfels = new R3Surfel [ nsurfels ];
    for (int i = 0; i < set->NPoints(); i++) {
      const R3SurfelPoint *point = set->Point(i);
      R3Point position = point->Position() - position_origin.Vector();
      R3Vector normal = point->Normal();
      R3Vector tangent = point->Tangent();
      R3Surfel *surfel = &surfels[i];
      surfel->SetPosition(position[0], position[1], position[2]);
      surfel->SetNormal(normal[0], normal[1], normal[2]);
      surfel->SetTangent(tangent[0], tangent[1], tangent[2]);
      surfel->SetRadius(0, point->Radius(0));
      surfel->SetRadius(1, point->Radius(1));
      surfel->SetTimestamp(point->Timestamp() - timestamp_origin);
      surfel->SetIdentifier(point->Identifier());
      surfel->SetColor(point->Color());
      surfel->SetFlags(point->Flags() & ~R3_SURFEL_MARKED_FLAG);
    }
  }
}



R3SurfelBlock::
R3SurfelBlock(const R3Surfel *surfels, int nsurfels,
  const R3Point& position_origin, RNScalar timestamp_origin)
  : surfels(NULL),
    nsurfels(nsurfels),
    position_origin(position_origin),
    bbox(FLT_MAX,FLT_MAX,FLT_MAX,-FLT_MAX,-FLT_MAX,-FLT_MAX),
    timestamp_origin(timestamp_origin),
    timestamp_range(FLT_MAX,-FLT_MAX),
    min_identifier(UINT_MAX),
    max_identifier(0),
    resolution(0),
    flags(0),
    data(NULL),
    database(NULL),
    database_index(-1),
    file_surfels_offset(0),
    file_surfels_count(0),
    file_read_count(0),
    node(NULL),
    opengl_id(0)
{
  // Copy surfels
  if (nsurfels > 0) {
    this->surfels = new R3Surfel [ nsurfels ];
    for (int i = 0; i < nsurfels; i++) {
      this->surfels[i] = surfels[i];
    }
  }
}



R3SurfelBlock::
R3SurfelBlock(const RNArray<const R3Surfel *>& array,
  const R3Point& position_origin, RNScalar timestamp_origin)
  : surfels(NULL),
    nsurfels(array.NEntries()),
    position_origin(position_origin),
    bbox(FLT_MAX,FLT_MAX,FLT_MAX,-FLT_MAX,-FLT_MAX,-FLT_MAX),
    timestamp_origin(timestamp_origin),
    timestamp_range(FLT_MAX,-FLT_MAX),
    min_identifier(UINT_MAX),
    max_identifier(0),
    resolution(0),
    flags(0),
    data(NULL),
    database(NULL),
    database_index(-1),
    file_surfels_offset(0),
    file_surfels_count(0),
    file_read_count(0),
    node(NULL),
    opengl_id(0)
{
  // Copy surfels
  if (nsurfels > 0) {
    surfels = new R3Surfel [ nsurfels ];
    for (int i = 0; i < nsurfels; i++) {
      this->surfels[i] = *(array[i]);
    }
  }
}



R3SurfelBlock::
R3SurfelBlock(const R3Point *points, int npoints)
  : surfels(NULL),
    nsurfels(npoints),
    position_origin(R3zero_point),
    bbox(FLT_MAX,FLT_MAX,FLT_MAX,-FLT_MAX,-FLT_MAX,-FLT_MAX),
    timestamp_origin(0),
    timestamp_range(FLT_MAX,-FLT_MAX),
    min_identifier(UINT_MAX),
    max_identifier(0),
    resolution(0),
    flags(0),
    data(NULL),
    database(NULL),
    database_index(-1),
    file_surfels_offset(0),
    file_surfels_count(0),
    file_read_count(0),
    node(NULL),
    opengl_id(0)
{
  // Check number of points
  if (npoints == 0) return;

  // Compute position_origin
  position_origin = R3Centroid(npoints, (R3Point *) points);

  // Compute bounding box
  for (int i = 0; i < npoints; i++) {
    bbox.Union(points[i]);
  }

  // Copy surfels
  surfels = new R3Surfel [ nsurfels ];
  for (int i = 0; i < nsurfels; i++) {
    R3Point point = points[i];
    point -= position_origin.Vector();
    R3Surfel& surfel = this->surfels[i];
    surfel.SetPosition(point.X(), point.Y(), point.Z());
    surfel.SetRadius(bbox.DiagonalLength() / npoints);
    surfel.SetColor(128, 128, 128);
  }
}



R3SurfelBlock::
R3SurfelBlock(const RNArray<R3Point *>& points)
  : surfels(NULL),
    nsurfels(points.NEntries()),
    position_origin(R3zero_point),
    bbox(FLT_MAX,FLT_MAX,FLT_MAX,-FLT_MAX,-FLT_MAX,-FLT_MAX),
    timestamp_origin(0),
    timestamp_range(FLT_MAX,-FLT_MAX),
    resolution(0),
    flags(0),
    data(NULL),
    database(NULL),
    database_index(-1),
    file_surfels_offset(0),
    file_surfels_count(0),
    file_read_count(0),
    node(NULL),
    opengl_id(0)
{
  // Check points
  if (points.IsEmpty()) return;
  
  // Compute position_origin
  position_origin = R3Centroid(points);

  // Compute bounding box
  for (int i = 0; i < points.NEntries(); i++) {
    bbox.Union(*(points[i]));
  }

  // Copy surfels
  surfels = new R3Surfel [ nsurfels ];
  for (int i = 0; i < nsurfels; i++) {
    R3Point point = *(points[i]);
    point -= position_origin.Vector();
    R3Surfel& surfel = this->surfels[i];
    surfel.SetPosition(point.X(), point.Y(), point.Z());
    surfel.SetRadius(bbox.DiagonalLength() / points.NEntries());
    surfel.SetColor(128, 128, 128);
  }
}



R3SurfelBlock::
~R3SurfelBlock(void)
{
  // Safety check
  assert(file_read_count == 0);

  // Remove from node
  if (node) node->RemoveBlock(this);

  // Remove from database
  if (database) database->RemoveBlock(this);

  // Delete surfels
  if (surfels) delete [] surfels;

#ifdef DRAW_WITH_DISPLAY_LIST
  // Delete opengl display lists
  if (opengl_id > 0) glDeleteLists(opengl_id, 2);
#endif

#ifdef DRAW_WITH_VBO
  glDeleteBuffers(1, &opengl_id);
#endif
}



////////////////////////////////////////////////////////////////////////
// ASSIGNMENT OPERATOR
////////////////////////////////////////////////////////////////////////

R3SurfelBlock& R3SurfelBlock::
operator=(const R3SurfelBlock& block)
{
  // Delete old surfels
  if (this->surfels) delete this->surfels;
  this->surfels = NULL;

  // Copy properties
  this->nsurfels = block.nsurfels;
  this->position_origin = block.position_origin;
  this->bbox = block.bbox;
  this->timestamp_origin = block.timestamp_origin;
  this->timestamp_range = block.timestamp_range;
  this->min_identifier = block.min_identifier;
  this->max_identifier = block.max_identifier;
  this->resolution = block.resolution;
  this->flags = block.flags & R3_SURFEL_BLOCK_PROPERTY_FLAGS;
  this->data = NULL;
  this->database = NULL;
  this->database_index = -1;
  this->file_surfels_offset = 0;
  this->file_surfels_count = 0;
  this->file_read_count = 0;
  this->node = NULL;
  this->opengl_id = 0;

  // Copy surfels
  if (this->nsurfels > 0) {
    this->surfels = new R3Surfel [ this->nsurfels ];
    for (int i = 0; i < this->nsurfels; i++) {
      this->surfels[i] = block.surfels[i];
    }
  }

  // Return this
  return *this;
}



////////////////////////////////////////////////////////////////////////
// PROPERTY FUNCTIONS
////////////////////////////////////////////////////////////////////////

const R3Box& R3SurfelBlock::
BBox(void) const
{
  // Update bbox
  if (!flags[R3_SURFEL_BLOCK_BBOX_UPTODATE_FLAG]) {
    R3SurfelBlock *block = (R3SurfelBlock *) this;
    block->UpdateBBox();
  }

  // Return bbox
  return bbox;
}



const RNInterval& R3SurfelBlock::
TimestampRange(void) const
{
  // Update timestamp range
  if (timestamp_range.IsEmpty()) {
    R3SurfelBlock *block = (R3SurfelBlock *) this;
    block->UpdateTimestampRange();
  }

  // Return timestamp range
  return timestamp_range;
}



unsigned int R3SurfelBlock::
MinIdentifier(void) const
{
  // Update min identifier
  if (max_identifier < min_identifier) {
    R3SurfelBlock *block = (R3SurfelBlock *) this;
    block->UpdateIdentifierRange();
  }

  // Return min identifier
  return min_identifier;
}



unsigned int R3SurfelBlock::
MaxIdentifier(void) const
{
  // Update max identifier
  if (max_identifier < min_identifier) {
    R3SurfelBlock *block = (R3SurfelBlock *) this;
    block->UpdateIdentifierRange();
  }

  // Return max identifier
  return max_identifier;
}



RNScalar R3SurfelBlock::
Resolution(void) const
{
  // Update resolution
  if (!flags[R3_SURFEL_BLOCK_RESOLUTION_UPTODATE_FLAG]) {
    R3SurfelBlock *block = (R3SurfelBlock *) this;
    block->UpdateResolution();
  }

  // Return resolution
  return resolution;
}



RNBoolean R3SurfelBlock::
HasActive(void) const
{
  // Update active flags
  if (!flags[R3_SURFEL_BLOCK_FLAGS_UPTODATE_FLAG]) {
    R3SurfelBlock *block = (R3SurfelBlock *) this;
    block->UpdateFlags();
  }

   // Return whether block has active 
  return flags[R3_SURFEL_BLOCK_HAS_ACTIVE_FLAG];
}



RNBoolean R3SurfelBlock::
HasNormals(void) const
{
  // Update normals flags
  if (!flags[R3_SURFEL_BLOCK_FLAGS_UPTODATE_FLAG]) {
    R3SurfelBlock *block = (R3SurfelBlock *) this;
    block->UpdateFlags();
  }

   // Return whether block has normals 
  return flags[R3_SURFEL_BLOCK_HAS_NORMALS_FLAG];
}



RNBoolean R3SurfelBlock::
HasTangents(void) const
{
  // Update tangents flags
  if (!flags[R3_SURFEL_BLOCK_FLAGS_UPTODATE_FLAG]) {
    R3SurfelBlock *block = (R3SurfelBlock *) this;
    block->UpdateFlags();
  }

   // Return whether block has tangents 
  return flags[R3_SURFEL_BLOCK_HAS_TANGENTS_FLAG];
}



RNBoolean R3SurfelBlock::
HasAerial(void) const
{
  // Update aerial/terrestrial flags
  if (!flags[R3_SURFEL_BLOCK_FLAGS_UPTODATE_FLAG]) {
    R3SurfelBlock *block = (R3SurfelBlock *) this;
    block->UpdateFlags();
  }

   // Return whether block has aerial scanner points
  return flags[R3_SURFEL_BLOCK_HAS_AERIAL_FLAG];
}



RNBoolean R3SurfelBlock::
HasTerrestrial(void) const
{
  // Update aerial/terrestrial flags
  if (!flags[R3_SURFEL_BLOCK_FLAGS_UPTODATE_FLAG]) {
    R3SurfelBlock *block = (R3SurfelBlock *) this;
    block->UpdateFlags();
  }

  // Return whether block has terrestrial scanner points
  return flags[R3_SURFEL_BLOCK_HAS_TERRESTRIAL_FLAG];
}



RNBoolean R3SurfelBlock::
IsDirty(void) const
{
  // Return whether block is dirty
  return flags[R3_SURFEL_BLOCK_DIRTY_FLAG];
}



void R3SurfelBlock::
SetDirty(RNBoolean dirty)
{
  // Set whether block is dirty
  if (dirty) flags.Add(R3_SURFEL_BLOCK_DIRTY_FLAG);
  else flags.Remove(R3_SURFEL_BLOCK_DIRTY_FLAG);
}



void R3SurfelBlock::
SetData(void *data) 
{
  // Set user data
  this->data = data;
}


  
////////////////////////////////////////////////////////////////////////
// PROPERTY MANIPULATION FUNCTIONS
////////////////////////////////////////////////////////////////////////

void R3SurfelBlock::
SetPositionOrigin(const R3Point& position)
{
  // Set position origin
  this->position_origin = position;
}



void R3SurfelBlock::
SetTimestampOrigin(RNScalar timestamp)
{
  // Set timestamp origin
  this->timestamp_origin = timestamp;
}



void R3SurfelBlock::
SetSurfelPosition(int surfel_index, const R3Point& position)
{
  // Check if surfels are resident
  if (!surfels) {
    RNFail("Unable to set surfel position for non-resident block\n");
    abort();
  }

  // Set surfel position (relative to position origin)
  float x = position[0] - position_origin[0];
  float y = position[1] - position_origin[1];
  float z = position[2] - position_origin[2];
  surfels[surfel_index].SetPosition(x, y, z);

  // Mark properties out of date
  flags.Remove(R3_SURFEL_BLOCK_BBOX_UPTODATE_FLAG);
  flags.Remove(R3_SURFEL_BLOCK_RESOLUTION_UPTODATE_FLAG);

  // Remember that block is dirty
  SetDirty();
}



void R3SurfelBlock::
SetSurfelNormal(int surfel_index, const R3Vector& normal)
{
  // Check if surfels are resident
  if (!surfels) {
    RNFail("Unable to set surfel position for non-resident block\n");
    abort();
  }

  // Set surfel normal
  surfels[surfel_index].SetNormal(normal.X(), normal.Y(), normal.Z());

  // Remember that block is dirty
  SetDirty();
}



void R3SurfelBlock::
SetSurfelTangent(int surfel_index, const R3Vector& tangent)
{
  // Check if surfels are resident
  if (!surfels) {
    RNFail("Unable to set surfel position for non-resident block\n");
    abort();
  }

  // Set surfel tangent
  surfels[surfel_index].SetTangent(tangent.X(), tangent.Y(), tangent.Z());

  // Remember that block is dirty
  SetDirty();
}



void R3SurfelBlock::
SetSurfelRadius(int surfel_index, RNLength radius)
{
  // Check if surfels are resident
  if (!surfels) {
    RNFail("Unable to set surfel position for non-resident block\n");
    abort();
  }

  // Set surfel normal
  surfels[surfel_index].SetRadius(radius);

  // Remember that block is dirty
  SetDirty();
}



void R3SurfelBlock::
SetSurfelRadius(int surfel_index, int axis, RNLength radius)
{
  // Check if surfels are resident
  if (!surfels) {
    RNFail("Unable to set surfel position for non-resident block\n");
    abort();
  }

  // Set surfel normal
  surfels[surfel_index].SetRadius(axis, radius);

  // Remember that block is dirty
  SetDirty();
}



void R3SurfelBlock::
SetSurfelColor(int surfel_index, const RNRgb& color)
{
  // Check if surfels are resident
  if (!surfels) {
    RNFail("Unable to set surfel position for non-resident block\n");
    abort();
  }

  // Set surfel color
  surfels[surfel_index].SetColor(color);

  // Remember that block is dirty
  SetDirty();
}



void R3SurfelBlock::
SetSurfelTimestamp(int surfel_index, RNScalar timestamp)
{
  // Check if surfels are resident
  if (!surfels) {
    RNFail("Unable to set surfel position for non-resident block\n");
    abort();
  }

  // Set surfel timestamp (relative to timestamp_origin)
  surfels[surfel_index].SetTimestamp(timestamp - timestamp_origin);

  // Remember that timestamp range is out of date
  timestamp_range.Empty();

  // Remember that block is dirty
  SetDirty();
}



void R3SurfelBlock::
SetSurfelIdentifier(int surfel_index, unsigned int identifier)
{
  // Check if surfels are resident
  if (!surfels) {
    RNFail("Unable to set surfel position for non-resident block\n");
    abort();
  }

  // Set surfel identifier
  surfels[surfel_index].SetIdentifier(identifier);

  // Remember that identifier range is out of date
  min_identifier = UINT_MAX;
  max_identifier = 0;

  // Remember that block is dirty
  SetDirty();
}



void R3SurfelBlock::
SetSurfelActive(int surfel_index, RNBoolean active)
{
  // Check if surfels are resident
  if (!surfels) {
    RNFail("Unable to set surfel position for non-resident block\n");
    abort();
  }

  // Set surfel active
  surfels[surfel_index].SetActive(active);

  // Remember that block is dirty
  SetDirty();
}



void R3SurfelBlock::
SetSurfelAerial(int surfel_index, RNBoolean aerial)
{
  // Check if surfels are resident
  if (!surfels) {
    RNFail("Unable to set surfel position for non-resident block\n");
    abort();
  }

  // Set surfel aerial
  surfels[surfel_index].SetAerial(aerial);

  // Remember that block is dirty
  SetDirty();
}



void R3SurfelBlock::
SetSurfelFlags(int surfel_index, unsigned char flags)
{
  // Check if surfels are resident
  if (!surfels) {
    RNFail("Unable to set surfel property for non-resident block\n");
    abort();
  }

  // Set surfel mark
  surfels[surfel_index].SetFlags(flags);

  // Remember that block is dirty
  SetDirty();
}



void R3SurfelBlock::
SetSurfelSilhouetteBoundary(int surfel_index, RNBoolean boundary)
{
  // Check if surfels are resident
  if (!surfels) {
    RNFail("Unable to set surfel property for non-resident block\n");
    abort();
  }

  // Set surfel mark
  surfels[surfel_index].SetSilhouetteBoundary(boundary);

  // Remember that block is dirty
  SetDirty();
}



void R3SurfelBlock::
SetSurfelShadowBoundary(int surfel_index, RNBoolean boundary)
{
  // Check if surfels are resident
  if (!surfels) {
    RNFail("Unable to set surfel property for non-resident block\n");
    abort();
  }

  // Set surfel mark
  surfels[surfel_index].SetShadowBoundary(boundary);

  // Remember that block is dirty
  SetDirty();
}



void R3SurfelBlock::
SetSurfelBorderBoundary(int surfel_index, RNBoolean boundary)
{
  // Check if surfels are resident
  if (!surfels) {
    RNFail("Unable to set surfel property for non-resident block\n");
    abort();
  }

  // Set surfel mark
  surfels[surfel_index].SetBorderBoundary(boundary);

  // Remember that block is dirty
  SetDirty();
}



void R3SurfelBlock::
SetSurfelMark(int surfel_index, RNBoolean mark)
{
  // Check if surfels are resident
  if (!surfels) {
    RNFail("Unable to set surfel property for non-resident block\n");
    abort();
  }

  // Set surfel mark
  surfels[surfel_index].SetMark(mark);
}



void R3SurfelBlock::
SetMarks(RNBoolean mark)
{
  // Check surfels
  if (!mark && !surfels) return;

  // Set mark for all surfels
  for (int i = 0; i < nsurfels; i++) {
    surfels[i].SetMark(mark);
  }
}



void R3SurfelBlock::
Transform(const R3Affine& transformation) 
{
  // Check if transformation is identity
  if (transformation.IsIdentity()) return;

  // Get scale factor
  RNScalar scale = transformation.ScaleFactor();
  if (RNIsEqual(scale, 1.0)) scale = 1.0;
  if (scale == 0) return;

  // Update resolution
  if (resolution > 0) resolution /= scale * scale;

  // Transform position_origin
  R3Point old_position_origin = position_origin;
  position_origin.Transform(transformation);

  // Read block
  if (database) database->ReadBlock(this);

  // Transform surfels 
  bbox = R3null_box;
  for (int i = 0; i < NSurfels(); i++) {
    R3Surfel *surfel = &surfels[i];
    R3Point position(surfel->X() + old_position_origin[0], surfel->Y() + old_position_origin[1], surfel->Z() + old_position_origin[2]);
    R3Vector normal(surfel->NX(), surfel->NY(), surfel->NZ());
    R3Vector tangent(surfel->TX(), surfel->TY(), surfel->TZ());
    position.Transform(transformation);
    normal.Transform(transformation);
    tangent.Transform(transformation);
    surfel->SetPosition(position.X() - position_origin.X(), position.Y() - position_origin.Y(), position.Z() - position_origin.Z());
    surfel->SetNormal(normal.X(), normal.Y(), normal.Z());
    surfel->SetTangent(tangent.X(), tangent.Y(), tangent.Z());
    surfel->SetRadius(0, scale * surfel->Radius(0));
    surfel->SetRadius(1, scale * surfel->Radius(1));
    bbox.Union(position);
  }

  // Remember that block is dirty
  SetDirty();

  // Release block
  if (database) database->ReleaseBlock(this);

  // Update database
  // ???
}



////////////////////////////////////////////////////////////////////////
// PROPERTY UPDATE FUNCTIONS
////////////////////////////////////////////////////////////////////////

void R3SurfelBlock::
UpdateProperties(void)
{
  // Read block
  if (database) database->ReadBlock(this);

  // Update properties
  double dummy = 0;
  dummy += BBox().Min().X();
  dummy += TimestampRange().Min();
  dummy += MinIdentifier();
  dummy += MaxIdentifier();
  dummy += Resolution();
  dummy += (HasAerial()) ? 1 : 2;
  if (dummy == 927612.21242) {
    printf("Amazing!\n");
  }

  // Update normals
  UpdateSurfelNormals();

  // Release block
  if (database) database->ReleaseBlock(this);
}



void R3SurfelBlock::
UpdateBBox(void)
{
  // Read block
  if (database) database->ReadBlock(this);

  // Update bounding box
  bbox = R3null_box;
  for (int i = 0; i < nsurfels; i++) {
    const float *p = surfels[i].PositionPtr();
    bbox.Union(R3Point(p[0], p[1], p[2]));
  }

  // Release block
  if (database) database->ReleaseBlock(this);

  // Translate bounding box by position_origin
  bbox.Translate(position_origin.Vector());

  // Mark bbox uptodate
  flags.Add(R3_SURFEL_BLOCK_BBOX_UPTODATE_FLAG);
}



void R3SurfelBlock::
UpdateTimestampRange(void)
{
  // Read block
  if (database) database->ReadBlock(this);

  // Update timestamp range
  timestamp_range = RNnull_interval;
  for (int i = 0; i < nsurfels; i++) {
    float t = surfels[i].Timestamp();
    timestamp_range.Union(t);
  }

  // Release block
  if (database) database->ReleaseBlock(this);

  // Translate bounding box by timestamp origin
  timestamp_range += timestamp_origin;
}



void R3SurfelBlock::
UpdateIdentifierRange(void)
{
  // Read block
  if (database) database->ReadBlock(this);

  // Update max identifier
  min_identifier = UINT_MAX;
  max_identifier = 0;
  for (int i = 0; i < nsurfels; i++) {
    unsigned int id = surfels[i].Identifier();
    if (id < min_identifier) min_identifier = id;
    if (id > max_identifier) max_identifier = id;
  }

  // Release block
  if (database) database->ReleaseBlock(this);
}



void R3SurfelBlock::
UpdateResolution(void)
{
  // Initialize resolution
  resolution = 0;

#if 0
  // Estimate resoluton based on surfel radii
  if (nsurfels > 0) {
    // Read block
    if (database) database->ReadBlock(this);

    // Sum surfel radii
    int nsamples = 1000;
    if (nsurfels < nsamples) nsamples = nsurfels;
    int step = nsurfels / nsamples;
    RNLength total_radius = 0;
    for (int i = 0; i < nsamples; i++) {
      R3Surfel *surfel = &surfels[i*step];
      RNLength radius = surfel->Radius();
      total_radius += radius;
    }

    // Resolution is samples / area
    if (RNIsZero(total_radius)) return;
    RNLength average_radius = total_radius / nsamples;
    RNScalar area = RN_PI * average_radius * average_radius;
    resolution = 1.0 / area;

    // Release block
    if (database) database->ReleaseBlock(this);
  }
#else
  // Estimate resolution based on number of surfels per cross-section of bounding box
  R3Box box = BBox();
  int dim = box.ShortestAxis();
  RNLength length1 = box.AxisLength((dim+1)%3);
  RNLength length2 = box.AxisLength((dim+2)%3);
  RNArea area = length1 * length2;
  resolution = (area > 0) ? sqrt(nsurfels / area) : 0;
#endif

  // Mark resolution uptodate
  flags.Add(R3_SURFEL_BLOCK_RESOLUTION_UPTODATE_FLAG);
}



void R3SurfelBlock::
UpdateFlags(void)
{
  // Read block
  if (database) database->ReadBlock(this);

  // Reset flags
  flags.Remove(R3_SURFEL_BLOCK_HAS_ACTIVE_FLAG);
  flags.Remove(R3_SURFEL_BLOCK_HAS_NORMALS_FLAG);
  flags.Remove(R3_SURFEL_BLOCK_HAS_AERIAL_FLAG);
  flags.Remove(R3_SURFEL_BLOCK_HAS_TERRESTRIAL_FLAG);

  // Update flags
  for (int i = 0; i < nsurfels; i++) {
    if (surfels[i].IsActive()) flags.Add(R3_SURFEL_BLOCK_HAS_ACTIVE_FLAG);
    if (surfels[i].HasNormal()) flags.Add(R3_SURFEL_BLOCK_HAS_NORMALS_FLAG);
    if (surfels[i].HasTangent()) flags.Add(R3_SURFEL_BLOCK_HAS_TANGENTS_FLAG);
    if (surfels[i].IsAerial()) flags.Add(R3_SURFEL_BLOCK_HAS_AERIAL_FLAG);
    else flags.Add(R3_SURFEL_BLOCK_HAS_TERRESTRIAL_FLAG);
  }

  // Release block
  if (database) database->ReleaseBlock(this);

  // Mark flags uptodate
  flags.Add(R3_SURFEL_BLOCK_FLAGS_UPTODATE_FLAG);
}



////////////////////////////////////////////////////////////////////////
// SURFEL UPDATE FUNCTIONS
////////////////////////////////////////////////////////////////////////

void R3SurfelBlock::
UpdateSurfelNormals(void) 
{
  // Update normals of all surfels in block
  R3SurfelPointSet pointset(this);
  pointset.UpdateNormals();
}



////////////////////////////////////////////////////////////////////////
// DRAW FUNCTIONS
////////////////////////////////////////////////////////////////////////

void R3SurfelBlock::
Draw(RNFlags flags) const
{
  // Get convenient variables
  int c = flags[R3_SURFEL_COLOR_DRAW_FLAG];
  int n = flags[R3_SURFEL_NORMAL_DRAW_FLAG];
  int id = flags[R3_SURFEL_IDENTIFIER_DRAW_FLAG];

  // Push translation to position_origin
  glPushMatrix();
  glTranslated(position_origin[0], position_origin[1], position_origin[2]);

#ifdef R3_SURFEL_DRAW_WITH_DISPLAY_LIST
  // Create a display list id to use to detect errors
  static GLuint error_id = 0;
  if (error_id == 0) error_id = glGenLists(1);

  // Create display list for block
  if (opengl_id == 0) {
    R3SurfelBlock *block = (R3SurfelBlock *) this;
    block->opengl_id = error_id;
    glGetError();
    GLuint id = glGenLists(2);
    if ((id > 0) && (glGetError() == GL_NO_ERROR)) {
      glNewList(id, GL_COMPILE);
#     if 1
        // Draw without color using arrays
        glEnableClientState(GL_VERTEX_ARRAY);
        glVertexPointer(3, GL_FLOAT, sizeof(R3Surfel), surfels);
        glDrawArrays(GL_POINTS, 0, NSurfels());
        glDisableClientState(GL_VERTEX_ARRAY);
#     else
        // Draw without color using glBegin and glEnd
        glBegin(GL_POINTS);
        for (int i = 0; i < NSurfels(); i++) {
          const R3Surfel& surfel = surfels[i];
          glVertex3fv(surfel.PositionPtr());
        }
        glEnd();
#     endif
      glEndList();
      if (glGetError() == GL_NO_ERROR) {
        glNewList(id+1, GL_COMPILE);
#       if 0
          // Draw with color using arrays
          glEnableClientState(GL_VERTEX_ARRAY);
          glEnableClientState(GL_COLOR_ARRAY);
          glVertexPointer(3, GL_FLOAT, sizeof(R3Surfel), surfels);
          glColorPointer(3, GL_UNSIGNED_BYTE, sizeof(R3Surfel), surfels[0].ColorPtr());
          glDrawArrays(GL_POINTS, 0, NSurfels());
          glDisableClientState(GL_VERTEX_ARRAY);
          glDisableClientState(GL_COLOR_ARRAY);
#       else
          // Draw with color using glBegin and glEnd
          glBegin(GL_POINTS);
          for (int i = 0; i < NSurfels(); i++) {
            const R3Surfel& surfel = surfels[i];
            glColor3ubv(surfel.ColorPtr());
            glVertex3fv(surfel.PositionPtr());
          }
          glEnd();
#       endif
        glEndList();
        if (glGetError() == GL_NO_ERROR) {
          block->opengl_id = id;
        }
      }
      if (opengl_id == error_id) {
        glDeleteLists(id, 2);
      }
    }
  }

  // Draw surfels
  if (opengl_id != error_id) {
    // Use display list if available
    glCallList((c) ? opengl_id + 1 : opengl_id);
  }
  else {
    // Draw surfels the slow way
    if (c) {
      // Draw surfels (with color)
      glBegin(GL_POINTS);
      for (int i = 0; i < NSurfels(); i++) {
        const R3Surfel& surfel = surfels[i];
        glColor3ubv(surfel.ColorPtr());
        glVertex3fv(surfel.PositionPtr());
      }
      glEnd();
    }
    else {
      // Draw surfels without color
      glBegin(GL_POINTS);
      for (int i = 0; i < NSurfels(); i++) {
        const R3Surfel& surfel = surfels[i];
        glVertex3fv(surfel.PositionPtr());
      }
      glEnd();
    }
  }
#endif

#ifdef R3_SURFEL_DRAW_WITH_VBO
  // Create a VBO id to use to detect errors
  static GLuint error_buffer_id = 0;
  if (error_buffer_id == 0) {
    glGenBuffers(1, &error_buffer_id);
  }

  // Create a VBO for block
  if (opengl_id == 0) {
    glGetError();
    GLuint buffer_id;
    opengl_id = error_buffer_id;
    glGenBuffers(1, &buffer_id);
    if (buffer_id > 0) {
      glBindBuffer(GL_ARRAY_BUFFER, buffer_id);
      glBufferData(GL_ARRAY_BUFFER, NSurfels() * sizeof(R3Surfel), surfels, GL_STATIC_DRAW);
      if (glGetError() == GL_NO_ERROR) {
        opengl_id = buffer_id;
      }
    }
  }

  // Draw surfels
  glEnableClientState(GL_VERTEX_ARRAY);
  if (c) glEnableClientState(GL_COLOR_ARRAY);
  if (opengl_id != error_buffer_id) {
    // Draw surfels using VBO arrays
    glBindBuffer(GL_ARRAY_BUFFER, opengl_id);
    glVertexPointer(3, GL_FLOAT, sizeof(R3Surfel), 0);
    if (c) glColorPointer(3, GL_UNSIGNED_BYTE, sizeof(R3Surfel), 12);
  }
  else {
    // Draw surfels using client-side arrays
    glVertexPointer(3, GL_FLOAT, sizeof(R3Surfel), surfels);
    if (c) glColorPointer(3, GL_UNSIGNED_BYTE, sizeof(R3Surfel), surfels[0].ColorPtr());
  }
  glDrawArrays(GL_POINTS, 0, NSurfels());
  glDisableClientState(GL_VERTEX_ARRAY);
  if (c) glDisableClientState(GL_COLOR_ARRAY);
#endif

#ifdef R3_SURFEL_DRAW_WITH_ARRAYS
  // Draw surfels using arrays
  glEnableClientState(GL_VERTEX_ARRAY);
  if (c) glEnableClientState(GL_COLOR_ARRAY);
  glVertexPointer(3, GL_FLOAT, sizeof(R3Surfel), surfels);
  if (c) glColorPointer(3, GL_UNSIGNED_BYTE, sizeof(R3Surfel), surfels[0].ColorPtr());
  glDrawArrays(GL_POINTS, 0, NSurfels());
  glDisableClientState(GL_VERTEX_ARRAY);
  if (c) glDisableClientState(GL_COLOR_ARRAY);
#endif

#ifdef R3_SURFEL_DRAW_WITH_GLBEGIN
  if (flags[R3_SURFEL_DISC_DRAW_FLAG] && HasNormals() && HasTangents()) {
    // Draw discs
    const int nsides = 6;
    glBegin(GL_TRIANGLES);
    for (int i = 0; i < NSurfels(); i++) {
      const R3Surfel& surfel = surfels[i];
      R3Vector normal(surfel.NX(), surfel.NY(), surfel.NZ());
      R3Vector tangent1(surfel.TX(), surfel.TY(), surfel.TZ());
      R3Vector tangent2 = normal % tangent1;
      double r1 = surfel.Radius(0);
      double r2 = surfel.Radius(1);
      if (r1 <= 0) r1 = 0.1;
      if (r2 <= 0) r2 = r1;
      if (c) glColor3ubv(surfel.ColorPtr());
      if (id) LoadUnsignedInt(surfel.Identifier());
      else if (n) R3LoadNormal(normal);
      R3Point p[nsides];
      for (int j = 0; j < nsides; j++) {
        double angle = RN_TWO_PI*j/nsides;
        p[j].Reset(surfel.PX(), surfel.PY(), surfel.PZ());
        p[j] += cos(angle) * r1 * tangent1;
        p[j] += sin(angle) * r2 * tangent2;
      }
      for (int j = 0; j < nsides; j++) {
        glVertex3fv(surfel.PositionPtr());
        R3LoadPoint(p[j]);
        R3LoadPoint(p[(j+1)%nsides]);
      }
    }
    glEnd();        
  }
  else {
    // Draw points
    glBegin(GL_POINTS);
    for (int i = 0; i < NSurfels(); i++) {
      const R3Surfel& surfel = surfels[i];
      if (c) glColor3ubv(surfel.ColorPtr());
      if (id) LoadUnsignedInt(surfel.Identifier());
      else if (n) glNormal3f(surfel.NX(), surfel.NY(), surfel.NZ());
      glVertex3fv(surfel.PositionPtr());
    }
    glEnd();
  }
#endif

  // Pop translation to position_origin
  glPopMatrix();
}



void R3SurfelBlock::
Print(FILE *fp, const char *prefix, const char *suffix) const
{
  // Check fp
  if (!fp) fp = stdout;

  // Print block
  if (prefix) fprintf(fp, "%s", prefix);
  fprintf(fp, "%d %d", DatabaseIndex(), NSurfels());
  if (suffix) fprintf(fp, "%s", suffix);
  fprintf(fp, "\n");
}



////////////////////////////////////////////////////////////////////////
// I/O FUNCTIONS
////////////////////////////////////////////////////////////////////////

int R3SurfelBlock::
ReadFile(const char *filename)
{
  // Parse input filename extension
  const char *extension;
  if (!(extension = strrchr(filename, '.'))) {
    printf("Filename %s has no extension (e.g., .ply)\n", filename);
    return 0;
  }

  // Read file of appropriate type
  if (!strncmp(extension, ".obj", 4)) {
    return ReadOBJFile(filename);
  }
  else if (!strncmp(extension, ".xyz", 4)) {
    return ReadXYZAsciiFile(filename);
  }
  else if (!strncmp(extension, ".bin", 4)) {
    return ReadXYZBinaryFile(filename);
  }
  else if (!strncmp(extension, ".blk", 4)) {
    return ReadBinaryFile(filename);
  }
  else if (!strncmp(extension, ".upc", 4)) {
    return ReadUPCFile(filename);
  }
  else { 
    RNFail("Unable to read file %s (unrecognized extension: %s)\n", filename, extension); 
    return 0; 
  }

  // Should never get here
  return 0;
}



int R3SurfelBlock::
WriteFile(const char *filename) const
{
  // Parse input filename extension
  const char *extension;
  if (!(extension = strrchr(filename, '.'))) {
    printf("Filename %s has no extension (e.g., .ply)\n", filename);
    return 0;
  }

  // Write file of appropriate type
  if (!strncmp(extension, ".xyz", 4)) {
    return WriteXYZAsciiFile(filename);
  }
  else if (!strncmp(extension, ".blk", 4)) {
    return WriteBinaryFile(filename);
  }
  else { 
    RNFail("Unable to write file %s (unrecognized extension: %s)\n", filename, extension); 
    return 0; 
  }

  // Should never get here
  return 0;
}



////////////////////////////////////////////////////////////////////////
// OBJ I/O FUNCTIONS
////////////////////////////////////////////////////////////////////////

int R3SurfelBlock::
ReadOBJFile(const char *filename)
{
  // Open file
  FILE *fp;
  if (!(fp = fopen(filename, "r"))) {
    RNFail("Unable to open file %s\n", filename);
    return 0;
  }

  // Read file
  if (!ReadOBJ(fp)) {
    RNFail("Unable to read OBJ file %s\n", filename);
    fclose(fp);
    return 0;
  }

  // Close file
  fclose(fp);

  // Return success
  return 1;
}



int R3SurfelBlock::
ReadOBJ(FILE *fp)
{
  // Get original file offset
  long int file_offset = RNFileTell(fp);

  // Count the number of surfels and compute centroid
  int count = 0;
  float cx = 0;
  float cy = 0;
  float cz = 0;
  char buffer[4096];
  while (fgets(buffer, 4096, fp)) {
    if (buffer[0] == 'v') {
      // Parse surfel data
      char keyword[64];
      float x, y, z;
      if (sscanf(buffer, "%s%f%f%f", keyword, &x, &y, &z) == (unsigned int) 4) {
        cx += x;
        cy += y;
        cz += z;
        count++;
      }
    }
  }

  // Check number of points
  if (count == 0) return 0;

  // Comput centroid of points
  cx /= count;
  cy /= count;
  cz /= count;

  // Set position_origin to centroid of points
  position_origin.Reset(cx, cy, cz);

  // Rewind file to original file offset
  RNFileSeek(fp, file_offset, RN_FILE_SEEK_SET);

  // Allocate surfels
  surfels = new R3Surfel [ count ];
  if (!surfels) {
    RNFail("Unable to allocate surfel block\n");
    return 0;
  }

  // Read surfels
  nsurfels = 0;
  while (fgets(buffer, 4096, fp)) {
    // Check number of surfels
    if (nsurfels >= count) break;

    // Check if point
    if (buffer[0] == 'v') {
      // Parse surfel data
      char keyword[64];
      float x, y, z;
      if (sscanf(buffer, "%s%f%f%f", keyword, &x, &y, &z) != (unsigned int) 4) {
        RNFail("Unable to read point %d out of %d into surfel block\n", nsurfels, count);
        delete [] surfels;
        surfels = NULL;
        nsurfels = 0;
        return 0;
      }

      // Assign surfel
      surfels[nsurfels].SetPosition(x - cx, y - cy, z - cz);
      nsurfels++;
    }
  }

  // Return success
  return 1;
}



////////////////////////////////////////////////////////////////////////
// XYZ I/O FUNCTIONS
////////////////////////////////////////////////////////////////////////

int R3SurfelBlock::
ReadXYZAsciiFile(const char *filename)
{
  // Open file
  FILE *fp;
  if (!(fp = fopen(filename, "r"))) {
    RNFail("Unable to open file %s\n", filename);
    return 0;
  }

  // Read file
  if (!ReadXYZAscii(fp)) {
    RNFail("Unable to read XYZ file %s\n", filename);
    fclose(fp);
    return 0;
  }

  // Close file
  fclose(fp);

  // Return success
  return 1;
}



int R3SurfelBlock::
WriteXYZAsciiFile(const char *filename) const
{
  // Open file
  FILE *fp;
  if (!(fp = fopen(filename, "w"))) {
    RNFail("Unable to open file %s\n", filename);
    return 0;
  }

  // Write file
  if (!WriteXYZAscii(fp)) {
    RNFail("Unable to write XYZ file %s\n", filename);
    fclose(fp);
    return 0;
  }

  // Close file
  fclose(fp);

  // Return success
  return 1;
}



int R3SurfelBlock::
ReadXYZAscii(FILE *fp)
{
  // Get original file offset
  long int file_offset = RNFileTell(fp);

  // Count the number of surfels
  int count = 0;
  char buffer[4096];
  while (fgets(buffer, 4096, fp)) {
    // Check if blank line
    char *bufferp = buffer;
    while (*bufferp && isspace(*bufferp)) bufferp++;
    if (!*bufferp) continue;
    count++;
  }

  // Rewind file to original file offset
  RNFileSeek(fp, file_offset, RN_FILE_SEEK_SET);

  // Allocate surfels
  surfels = new R3Surfel [ count ];
  if (!surfels) {
    RNFail("Unable to allocate surfel block\n");
    return 0;
  }

  // Read surfels
  nsurfels = 0;
  while (fgets(buffer, 4096, fp)) {
    // Check if blank line
    char *bufferp = buffer;
    while (*bufferp && isspace(*bufferp)) bufferp++;
    if (!*bufferp) continue;

    // Check number of surfels
    if (nsurfels >= count) break;

    // Parse surfel data
    float x, y, z;
    unsigned int r, g, b;
    if (sscanf(buffer, "%f%f%f%u%u%u", &x, &y, &z, &r, &g, &b) != (unsigned int) 6) {
      r = 255; g = 0; b = 0;
      if (sscanf(buffer, "%f%f%f", &x, &y, &z) != (unsigned int) 3) {
        RNFail("Unable to read point %d out of %d into surfel block\n", nsurfels, count);
        delete [] surfels;
        surfels = NULL;
        nsurfels = 0;
        return 0;
      }
    }

    // Assign surfel
    surfels[nsurfels].SetPosition(x, y, z);
    surfels[nsurfels].SetColor(r, g, b);
    nsurfels++;
  }

  // Return success
  return 1;
}



int R3SurfelBlock::
WriteXYZAscii(FILE *fp) const
{
  // Write surfels
  for (int i = 0; i < NSurfels(); i++) {
    const R3Surfel *surfel = Surfel(i);
    const float *position = surfel->PositionPtr();
    const unsigned char *color = surfel->ColorPtr();
    fprintf(fp, "%g %g %g ", position[0], position[1], position[2]);
    fprintf(fp, "%u %u %u ", color[0], color[1], color[2]);
    fprintf(fp, "0 0\n");
  }

  // Return success
  return 1;
}



////////////////////////////////////////////////////////////////////////
// FLOAT I/O FUNCTIONS
////////////////////////////////////////////////////////////////////////

int R3SurfelBlock::
ReadXYZBinaryFile(const char *filename)
{
  // Open file
  FILE *fp;
  if (!(fp = fopen(filename, "rb"))) {
    RNFail("Unable to open file %s\n", filename);
    return 0;
  }

  // Read file
  if (!ReadXYZBinary(fp)) {
    RNFail("Unable to read XYZ binary file %s\n", filename);
    fclose(fp);
    return 0;
  }

  // Close file
  fclose(fp);

  // Return success
  return 1;
}



int R3SurfelBlock::
ReadXYZBinary(FILE *fp)
{
  // Determine the number of surfels
  long int start_file_offset = RNFileTell(fp);
  RNFileSeek(fp, 0, RN_FILE_SEEK_END);
  long int end_file_offset = RNFileTell(fp);
  RNFileSeek(fp, start_file_offset, RN_FILE_SEEK_SET);
  long int file_size = end_file_offset - start_file_offset;
  nsurfels = file_size / (4*sizeof(RNScalar32));
  if (nsurfels == 0) return 0;
  
  // Allocate surfels
  surfels = new R3Surfel [ nsurfels ];
  if (!surfels) {
    RNFail("Unable to allocate surfel block\n");
    return 0;
  }

  // Read surfels
  position_origin = R3zero_point;
  for (int i = 0; i < nsurfels; i++) {
    // Read data
    float xyzr[4];
    if (fread(xyzr, sizeof(RNScalar32), 4, fp) != (size_t) 4) {
      RNFail("Unable to read point %d\n", i);
      delete [] surfels;
      surfels = NULL;
      nsurfels = 0;
      return 0;
    }

    // Compute color
    float c = 255 * xyzr[3];

    // Assign surfel position
    surfels[i].SetPosition(xyzr);
    surfels[i].SetColor(c, c, c);

    // Update position_origin
    position_origin[0] += xyzr[0];
    position_origin[1] += xyzr[1];
    position_origin[2] += xyzr[2];
  }

  // Put position_origin at centroid
  position_origin /= nsurfels;

  // Update surfels to be relative to position_origin
  for (int i = 0; i < nsurfels; i++) {
    const float *xyz = surfels[i].PositionPtr();
    surfels[i].SetPosition(xyz[0] - position_origin[0], xyz[1] - position_origin[1], xyz[2] - position_origin[2]);
  }  
  
  // Return success
  return 1;
}



////////////////////////////////////////////////////////////////////////
// BINARY I/O FUNCTIONS
////////////////////////////////////////////////////////////////////////

int R3SurfelBlock::
ReadBinaryFile(const char *filename)
{
  // Open file
  FILE *fp;
  if (!(fp = fopen(filename, "rb"))) {
    RNFail("Unable to open file %s\n", filename);
    return 0;
  }

  // Read file
  if (!ReadBinary(fp)) {
    RNFail("Unable to read surfel file %s\n", filename);
    fclose(fp);
    return 0;
  }

  // Close file
  fclose(fp);

  // Return success
  return 1;
}



int R3SurfelBlock::
WriteBinaryFile(const char *filename) const
{
  // Open file
  FILE *fp;
  if (!(fp = fopen(filename, "wb"))) {
    RNFail("Unable to open file %s\n", filename);
    return 0;
  }

  // Write file
  if (!WriteBinary(fp)) {
    RNFail("Unable to write surfel file %s\n", filename);
    fclose(fp);
    return 0;
  }

  // Close file
  fclose(fp);

  // Return success
  return 1;
}



int R3SurfelBlock::
ReadBinary(FILE *fp)
{
  // Read number of surfels
  if (fread(&nsurfels, sizeof(int), 1, fp) != (size_t) 1) {
    RNFail("Unable to read number of surfels\n");
    return 0;
  }

  // Read bounding box
  if (fread(&bbox, sizeof(R3Box), 1, fp) != (size_t) 1) {
    RNFail("Unable to read bounding box of surfel block\n");
    return 0;
  }

  // Read resolution
  if (fread(&resolution, sizeof(RNLength), 1, fp) != (size_t) 1) {
    RNFail("Unable to read resolution of surfel block\n");
    return 0;
  }

  // Read flags
  if (fread(&flags, sizeof(RNFlags), 1, fp) != (size_t) 1) {
    RNFail("Unable to read flags of surfel block\n");
    return 0;
  }

  // Allocate memory for surfels
  surfels = new R3Surfel [ nsurfels ];
  if (!surfels) {
    RNFail("Unable to allocate surfel block\n");
    return 0;
  }

  // Read surfels
  int count = 0;
  while (count < nsurfels) {
    int status = fread(&surfels[count], sizeof(R3Surfel), nsurfels - count, fp);
    if (status <= 0) {
      RNFail("Unable to read surfel block\n");
      delete [] surfels;
      surfels = NULL;
      nsurfels = 0;
      return 0;
    }
    count += status;
  }

  // Return success
  return 1;
}



int R3SurfelBlock::
WriteBinary(FILE *fp) const
{
  // Write number of surfels
  int nsurfels = NSurfels();
  if (fwrite(&nsurfels, sizeof(int), 1, fp) != (size_t) 1) {
    RNFail("Unable to write number of surfels\n");
    return 0;
  }

  // Write bounding box
  if (fwrite(&bbox, sizeof(R3Box), 1, fp) != (size_t) 1) {
    RNFail("Unable to write bounding box of surfel block\n");
    return 0;
  }

  // Write resolution
  if (fwrite(&resolution, sizeof(R3Box), 1, fp) != (size_t) 1) {
    RNFail("Unable to write resolution of surfel block\n");
    return 0;
  }

  // Write flags
  RNFlags tmp = flags; tmp.Remove(R3_SURFEL_BLOCK_DATABASE_FLAGS);
  if (fwrite(&tmp, sizeof(RNFlags), 1, fp) != (size_t) 1) {
    RNFail("Unable to write flags of surfel block\n");
    return 0;
  }

  // Write surfels
  for (int i = 0; i < NSurfels(); i++) {
    if (fwrite(&surfels[i], sizeof(R3Surfel), 1, fp) != (size_t) 1) {
      return 0;
    }
  }

  // Return success
  return 1;
}



////////////////////////////////////////////////////////////////////////
// UPC I/O FUNCTIONS
////////////////////////////////////////////////////////////////////////

struct UPCHeader {
  char signature[4]; // Signature of the file format (always UPCf)
  RNUChar8 versionMajor; // The major version number for the file
  RNUChar8 versionMinor; // The minor version number for the file
  RNUInt16 headerSize; // Size of the header block
  RNInt64 numOfPts; // The number of points within the file
  RNScalar64 xScale; // The scale used in the x-coord
  RNScalar64 yScale; // The scale used in the y-coord
  RNScalar64 zScale; // The scale used in the z-coord
  RNScalar64 xOffset; // The offset used in the x-coord
  RNScalar64 yOffset; // The offset used in the y-coord
  RNScalar64 zOffset; // The offset used in the z-coord
};

struct UPCPoint {
  RNScalar64 gpsTimeOfWeek; // The GPS time of week of the point
  RNUChar8 sensorNum[2]; // The laser sensor number used for the point
  RNUChar8 julianDay[3]; // The day the point was collected
  RNUChar8 flightLine[3]; // The flight line number of the point
  RNInt32 x; // The recorded x-coord of the point
  RNInt32 y; // The recorded y-coord of the point
  RNInt32 z; // The recorded z-coord of the point
  RNUChar8 intensity; // The intensity of the point
  RNUInt16 red; // The red component of the point
  RNUInt16 green; // The green component of the point
  RNUInt16 blue; // The blue component of the point
  RNUChar8 returnNum; // The return number of the point
};



int R3SurfelBlock::
ReadUPCFile(const char *filename)
{
  // Open file
  FILE *fp;
  if (!(fp = fopen(filename, "rb"))) {
    RNFail("Unable to open file %s\n", filename);
    return 0;
  }

  // Read file
  if (!ReadUPC(fp)) {
    RNFail("Unable to read surfel file %s\n", filename);
    return 0;
  }

  // Close file
  fclose(fp);

  // Return success
  return 1;
}



static int 
ReadUPCPreamble(FILE *fp)
{
  // Read preamble (some upc files have junk at front!)
  unsigned long offset = 0;
  char signature[8] = { '\0' };
  while (TRUE) {
    // Read/check signature
    if (fread(signature, sizeof(char), 4, fp) != (size_t) 4) return 0;
    if (!strcmp(signature, "UPCf")) break;
    else offset += 4;
  }
   
  // Seek to start of header
  RNFileSeek(fp, offset, RN_FILE_SEEK_SET);

  // Return success
  return 1;
}



static int 
ReadUPCHeader(FILE *fp, UPCHeader& header)
{
  // Read preamble
  if (!ReadUPCPreamble(fp)) {
    RNFail("Unable to read signature of UPC file\n");
    return 0;
  }

  // Read signature
  if (fread(&header.versionMajor, 1, 4, fp) != (size_t) 4) {
    RNFail("Unable to read UPC header signature\n");
    return 0;
  }

  // Read versionMajor
  if (fread(&header.versionMajor, 1, 1, fp) != (size_t) 1) {
    RNFail("Unable to read UPC header versionMajor\n");
    return 0;
  }

  // Read versionMinor
  if (fread(&header.versionMinor, 1, 1, fp) != (size_t) 1) {
    RNFail("Unable to read UPC header versionMinor\n");
    return 0;
  }

  // Read headerSize
  if (fread(&header.headerSize, 2, 1, fp) != (size_t) 1) {
    RNFail("Unable to read UPC header size\n");
    return 0;
  }

  // Read number of points
  if (fread(&header.numOfPts, 8, 1, fp) != (size_t) 1) {
    RNFail("Unable to read number of points in header of UPC file\n");
    return 0;
  }

  // Read scales and offsets
  if (fread(&header.xScale, 8, 6, fp) != (size_t) 6) {
    RNFail("Unable to read scales and offsets in header of UPC file\n");
    return 0;
  }

  // Return success
  return 1;
}



int R3SurfelBlock::
ReadUPC(FILE *fp)
{
  // Read header
  UPCHeader header;
  if (!ReadUPCHeader(fp, header)) return 0;

  // Compute scale and offset
  double xoffset = header.xOffset;
  double yoffset = header.yOffset;
  double zoffset = header.zOffset;
  double xscale = header.xScale;
  double yscale = header.yScale;
  double zscale = header.zScale;
 
  // XXX THIS IS A HACK FOR OTTAWA XXX
  xoffset -= 444965;
  yoffset -= 5029450;

  // Set position_origin
  position_origin = R3Point(xoffset, yoffset, zoffset);

  // Allocate memory for surfels
  int target_count = (int) header.numOfPts;
  surfels = new R3Surfel [ target_count ];
  if (!surfels) {
    RNFail("Unable to allocate surfel block\n");
    return 0;
  }

  // Allocate memory for UPC points
  const int upc_point_size = 36;
  char *upc_points = new char [ target_count * upc_point_size ];

  // Read upc points
  int upc_count = 0;
  while (upc_count < target_count) {
    int status = fread(&upc_points[upc_count], upc_point_size, target_count - upc_count, fp);
    if (status <= 0) break;
    upc_count += status;
  }

  // Assign surfels
  nsurfels = 0;
  for (int i = 0; i < upc_count; i++) {
    char *upc_point = &upc_points[i * upc_point_size];
    double x = *((RNInt32 *) &upc_point[16]) * xscale;
    double y = *((RNInt32 *) &upc_point[20]) * yscale;
    double z = *((RNInt32 *) &upc_point[24]) * zscale;
    RNUInt16 r = *((RNUInt16 *) &upc_point[29]);
    RNUInt16 g = *((RNUInt16 *) &upc_point[31]);
    RNUInt16 b = *((RNUInt16 *) &upc_point[33]);
    if (r > 255) r = 255;
    if (g > 255) g = 255;
    if (b > 255) b = 255;
    surfels[nsurfels].SetPosition(x, y, z);
    surfels[nsurfels].SetColor(r, g, b);
    surfels[nsurfels].SetAerial(upc_point[9] == '0');
    nsurfels++;
  }

  // Delete memory for UPC points
  delete [] upc_points;

  // Return success
  return 1;
}



////////////////////////////////////////////////////////////////////////
// INTERNAL FUNCTIONS (DO NOT USE)
////////////////////////////////////////////////////////////////////////

void R3SurfelBlock::
ResetSurfels(int nsurfels)
{
  // Delete old surfels
  if (surfels) delete [] surfels;

  // Reset everything
  this->surfels = NULL;
  this->nsurfels = 0;
  this->position_origin = R3zero_point;
  this->bbox = R3null_box;
  this->timestamp_origin = 0;
  this->timestamp_range.Reset(FLT_MAX,-FLT_MAX);
  this->min_identifier = UINT_MAX;
  this->max_identifier = 0;
  this->resolution = 0;
  this->flags = 0;
  this->data = NULL;
  this->database = NULL;
  this->database_index = -1;
  this->file_surfels_offset = 0;
  this->file_surfels_count = 0;
  this->file_read_count = 0;
  this->node = NULL;
  this->opengl_id = 0;

  // Allocate new surfels
  if (nsurfels > 0) {
    this->nsurfels = nsurfels;
    this->surfels = new R3Surfel [ nsurfels ];
  }
}



////////////////////////////////////////////////////////////////////////
// UPDATE FUNCTIONS
////////////////////////////////////////////////////////////////////////

void R3SurfelBlock::
UpdateAfterInsert(R3SurfelNode *node)
{
  // Update node
  this->node = node;
}



void R3SurfelBlock::
UpdateBeforeRemove(R3SurfelNode *node)
{
  // Update node
  this->node = NULL;
}



} // namespace gaps
