/* Source file for the R3 surfel constraint class */



////////////////////////////////////////////////////////////////////////
// INCLUDE FILES
////////////////////////////////////////////////////////////////////////

#include "R3Surfels.h"



////////////////////////////////////////////////////////////////////////
// Namespace
////////////////////////////////////////////////////////////////////////

namespace gaps {



////////////////////////////////////////////////////////////////////////
// BASE CONSTRAINT FUNCTIONS
////////////////////////////////////////////////////////////////////////

R3SurfelConstraint::
R3SurfelConstraint(void)
{
}



R3SurfelConstraint::
~R3SurfelConstraint(void)
{
}



int R3SurfelConstraint::
Check(const R3SurfelObject *object) const
{
  // Return PASS if all surfels in object satisfy constraint
  // Return MAYBE if some surfels in object may satisfy constraint
  // Return FAIL if no surfels in object could possibly satisfy constraint
  return Check(object->BBox());
}



int R3SurfelConstraint::
Check(const R3SurfelNode *node) const
{
  // Return PASS if all surfels in node satisfy constraint
  // Return MAYBE if some surfels in node may satisfy constraint
  // Return FAIL if no surfels in node could possibly satisfy constraint
  return Check(node->BBox());
}



int R3SurfelConstraint::
Check(const R3SurfelBlock *block) const
{
  // Return PASS if all surfels in block satisfy constraint
  // Return MAYBE if some surfels in block may satisfy constraint
  // Return FAIL if no surfels in block could possibly satisfy constraint
  return Check(block->BBox());
}



int R3SurfelConstraint::
Check(const R3SurfelBlock *block, const R3Surfel *surfel) const
{
  // Return PASS if surfel satisfies constraint
  // Return FAIL if surfel does not satisfy constraint
  return Check(block->PositionOrigin() + R3Point(surfel->X(), surfel->Y(), surfel->Z()));
}



int R3SurfelConstraint::
Check(const R3Box& box) const
{
  // Return PASS if all surfels in box satisfy constraint
  // Return MAYBE if some surfels in box may satisfy constraint
  // Return FAIL if no surfels in box could possibly satisfy constraint
  return R3_SURFEL_CONSTRAINT_MAYBE;
}



int R3SurfelConstraint::
Check(const R3Point& point) const
{
  // Return PASS if all surfels at this position satisfy constraint
  // Return MAYBE if some surfels at this position may satisfy constraint
  // Return FAIL if no surfels at this position could possibly satisfy constraint
  return R3_SURFEL_CONSTRAINT_PASS;
}



////////////////////////////////////////////////////////////////////////
// TIMESTAMP CONSTRAINT FUNCTIONS
////////////////////////////////////////////////////////////////////////

R3SurfelTimestampConstraint::
R3SurfelTimestampConstraint(const RNInterval& timestamp_range)
  : timestamp_range(timestamp_range)
{
}


  
R3SurfelTimestampConstraint::
R3SurfelTimestampConstraint(double min_timestamp, double max_timestamp)
  : timestamp_range(min_timestamp, max_timestamp)
{
}



int R3SurfelTimestampConstraint::
Check(const R3SurfelObject *object) const
{
  // Return whether object intersects timestamp range
  if (timestamp_range.IsEmpty()) return R3_SURFEL_CONSTRAINT_PASS;
  if (object->TimestampRange().IsEmpty()) return R3_SURFEL_CONSTRAINT_FAIL;
  if (timestamp_range.Intersects(object->TimestampRange())) return R3_SURFEL_CONSTRAINT_PASS;
  return R3_SURFEL_CONSTRAINT_FAIL;
}


  
int R3SurfelTimestampConstraint::
Check(const R3SurfelNode *node) const
{
  // Return whether node intersects timestamp range
  if (timestamp_range.IsEmpty()) return R3_SURFEL_CONSTRAINT_PASS;
  if (node->TimestampRange().IsEmpty()) return R3_SURFEL_CONSTRAINT_FAIL;
  if (timestamp_range.Intersects(node->TimestampRange())) return R3_SURFEL_CONSTRAINT_PASS;
  return R3_SURFEL_CONSTRAINT_FAIL;
}



int R3SurfelTimestampConstraint::
Check(const R3SurfelBlock *block) const
{
  // Return whether block intersects timestamp range
  if (timestamp_range.IsEmpty()) return R3_SURFEL_CONSTRAINT_PASS;
  if (block->TimestampRange().IsEmpty()) return R3_SURFEL_CONSTRAINT_FAIL;
  if (timestamp_range.Intersects(block->TimestampRange())) return R3_SURFEL_CONSTRAINT_PASS;
  return R3_SURFEL_CONSTRAINT_FAIL;
}



int R3SurfelTimestampConstraint::
Check(const R3SurfelBlock *block, const R3Surfel *surfel) const
{
  // Return whether block intersects timestamp range
  int surfel_index = block->SurfelIndex(surfel);
  assert((surfel_index >= 0) && (surfel_index < block->NSurfels()));
  if (timestamp_range.IsEmpty()) return R3_SURFEL_CONSTRAINT_PASS;
  if (timestamp_range.Contains(block->SurfelTimestamp(surfel_index))) return R3_SURFEL_CONSTRAINT_PASS;
  return R3_SURFEL_CONSTRAINT_FAIL;
}



////////////////////////////////////////////////////////////////////////
// COORDINATE CONSTRAINT FUNCTIONS
////////////////////////////////////////////////////////////////////////

R3SurfelCoordinateConstraint::
R3SurfelCoordinateConstraint(RNDimension dimension, const RNInterval& interval)
  : dimension(dimension), interval(interval)
{
}



int R3SurfelCoordinateConstraint::
Check(const R3Box& box) const
{
  // Return whether any point in box can satisfy constraint
  RNInterval box_interval(box.Min()[dimension], box.Max()[dimension]);
  if (interval.Contains(box_interval)) return R3_SURFEL_CONSTRAINT_PASS;
  else if (interval.Intersects(box_interval)) return R3_SURFEL_CONSTRAINT_MAYBE;
  return R3_SURFEL_CONSTRAINT_FAIL;
}



int R3SurfelCoordinateConstraint::
Check(const R3Point& point) const
{
  // Return whether coordinate value is set
  if (interval.Contains(point.Coord(dimension))) return R3_SURFEL_CONSTRAINT_PASS;
  return R3_SURFEL_CONSTRAINT_FAIL;
}



////////////////////////////////////////////////////////////////////////
// NORMAL CONSTRAINT FUNCTIONS
////////////////////////////////////////////////////////////////////////

R3SurfelNormalConstraint::
R3SurfelNormalConstraint(const R3Vector& direction, RNAngle max_angle)
{
  // Remember reference vector
  this->direction[0] = direction.X();
  this->direction[1] = direction.Y();
  this->direction[2] = direction.Z();

  // Compute min dot product
  this->min_dot = cos(max_angle);
}



int R3SurfelNormalConstraint::
Check(const R3SurfelBlock *block, const R3Surfel *surfel) const
{
  // Return whether surfel normal is within max_angle of reference
  RNScalar dot = 0.0;
  dot += surfel->NX() * direction[0];
  dot += surfel->NY() * direction[1];
  dot += surfel->NZ() * direction[2];
  if (dot < min_dot) return R3_SURFEL_CONSTRAINT_FAIL;
  return R3_SURFEL_CONSTRAINT_PASS;
}



////////////////////////////////////////////////////////////////////////
// BOX CONSTRAINT FUNCTIONS
////////////////////////////////////////////////////////////////////////

R3SurfelBoxConstraint::
R3SurfelBoxConstraint(const R3Box& box)
  : box(box)
{
}



int R3SurfelBoxConstraint::
Check(const R3Box& query) const
{
  // Return whether box contains/intersects query
  if (R3Contains(box, query)) return R3_SURFEL_CONSTRAINT_PASS;
  else if (R3Intersects(box, query)) return R3_SURFEL_CONSTRAINT_MAYBE;
  return R3_SURFEL_CONSTRAINT_FAIL;
}



int R3SurfelBoxConstraint::
Check(const R3Point& point) const
{
  // Return whether box contains point
  if (R3Contains(box, point)) return R3_SURFEL_CONSTRAINT_PASS;
  return R3_SURFEL_CONSTRAINT_FAIL;
}



////////////////////////////////////////////////////////////////////////
// CYLINDER CONSTRAINT FUNCTIONS
////////////////////////////////////////////////////////////////////////

R3SurfelCylinderConstraint::
R3SurfelCylinderConstraint(const R3Point& center, RNLength radius, RNCoord zmin, RNCoord zmax)
  : center(center),
    radius_squared(radius*radius),
    zmin(zmin),
    zmax(zmax)
{
}



int R3SurfelCylinderConstraint::
Check(const R3Box& box) const
{
  // Should check whether cylinder contains box
  // ???

  // Return whether cylinder intersects box
  RNScalar dx = 0, dy = 0;
  if (center.X() < box.XMin()) dx = box.XMin() - center.X();
  else if (center.X() > box.XMax()) dx = center.X() - box.XMax();
  if (center.Y() < box.YMin()) dy = box.YMin() - center.Y();
  else if (center.Y() > box.YMax()) dy = center.Y() - box.YMax();
  RNScalar d_squared = dx*dx + dy*dy;
  if (d_squared > radius_squared) return R3_SURFEL_CONSTRAINT_FAIL;
  if (box.ZMax() < zmin) return R3_SURFEL_CONSTRAINT_FAIL;
  if (box.ZMin() > zmax) return R3_SURFEL_CONSTRAINT_FAIL;
  return R3_SURFEL_CONSTRAINT_MAYBE;
}



int R3SurfelCylinderConstraint::
Check(const R3Point& point) const
{
  // Return whether cylinder contains point
  RNScalar dx = point.X() - center.X();
  RNScalar dy = point.Y() - center.Y();
  RNScalar d_squared = dx*dx + dy*dy;
  if (d_squared > radius_squared) return R3_SURFEL_CONSTRAINT_FAIL;
  if (point.Z() < zmin) return R3_SURFEL_CONSTRAINT_FAIL;
  if (point.Z() > zmax) return R3_SURFEL_CONSTRAINT_FAIL;
  return R3_SURFEL_CONSTRAINT_PASS;
}



////////////////////////////////////////////////////////////////////////
// SPHERE CONSTRAINT FUNCTIONS
////////////////////////////////////////////////////////////////////////

R3SurfelSphereConstraint::
R3SurfelSphereConstraint(const R3Sphere& sphere)
  : sphere(sphere)
{
}



int R3SurfelSphereConstraint::
Check(const R3Box& box) const
{
  // Return whether sphere contains/intersects box
  if (R3Contains(sphere, box)) return R3_SURFEL_CONSTRAINT_PASS;
  else if (R3Intersects(sphere, box)) return R3_SURFEL_CONSTRAINT_MAYBE;
  return R3_SURFEL_CONSTRAINT_FAIL;
}



int R3SurfelSphereConstraint::
Check(const R3Point& point) const
{
  // Return whether sphere contains point
  if (R3Contains(sphere, point)) return R3_SURFEL_CONSTRAINT_PASS;
  return R3_SURFEL_CONSTRAINT_FAIL;
}



////////////////////////////////////////////////////////////////////////
// HALFSPACE CONSTRAINT FUNCTIONS
////////////////////////////////////////////////////////////////////////

R3SurfelHalfspaceConstraint::
R3SurfelHalfspaceConstraint(const R3Halfspace& halfspace)
  : halfspace(halfspace)
{
}



int R3SurfelHalfspaceConstraint::
Check(const R3Box& box) const
{
  // Return whether halfspace contains/intersects box
  if (R3Contains(halfspace, box)) return R3_SURFEL_CONSTRAINT_PASS;
  else if (R3Intersects(halfspace, box)) return R3_SURFEL_CONSTRAINT_MAYBE;
  return R3_SURFEL_CONSTRAINT_FAIL;
}



int R3SurfelHalfspaceConstraint::
Check(const R3Point& point) const
{
  // Return whether halfspace contains point
  if (R3Contains(halfspace, point)) return R3_SURFEL_CONSTRAINT_PASS;
  return R3_SURFEL_CONSTRAINT_FAIL;
}



////////////////////////////////////////////////////////////////////////
// LINE CONSTRAINT FUNCTIONS
////////////////////////////////////////////////////////////////////////

R3SurfelLineConstraint::
R3SurfelLineConstraint(const R3Line& line, RNLength tolerance)
  : line(line),
    tolerance(tolerance)
{
}



int R3SurfelLineConstraint::
Check(const R3Box& box) const
{
  // Return whether contains/intersects box
  RNScalar d = R3Distance(line, box);
  if (d <= tolerance) return R3_SURFEL_CONSTRAINT_MAYBE;
  return R3_SURFEL_CONSTRAINT_FAIL;
}



int R3SurfelLineConstraint::
Check(const R3Point& point) const
{
  // Return whether line contains point
  RNScalar d = R3Distance(line, point);
  if (d <= tolerance) return R3_SURFEL_CONSTRAINT_PASS;
  return R3_SURFEL_CONSTRAINT_FAIL;
}



////////////////////////////////////////////////////////////////////////
// PLANE CONSTRAINT FUNCTIONS
////////////////////////////////////////////////////////////////////////

R3SurfelPlaneConstraint::
R3SurfelPlaneConstraint(const R3Plane& plane, RNBoolean below, RNBoolean on, RNBoolean above, RNLength tolerance)
  : plane(plane),
    below(below),
    on(on),
    above(above),
    tolerance(tolerance)
{
}



int R3SurfelPlaneConstraint::
Check(const R3Box& box) const
{
  // Returncontains/intersects box
  RNScalar d = R3SignedDistance(plane, box);
  if (fabs(d) <= tolerance) return R3_SURFEL_CONSTRAINT_MAYBE;
  if (above && (d > tolerance)) return R3_SURFEL_CONSTRAINT_PASS;
  if (below && (d < tolerance)) return R3_SURFEL_CONSTRAINT_PASS;
  return R3_SURFEL_CONSTRAINT_FAIL;
}



int R3SurfelPlaneConstraint::
Check(const R3Point& point) const
{
  // Return whether plane contains point
  RNScalar d = R3SignedDistance(plane, point);
  if (above && (d > tolerance)) return R3_SURFEL_CONSTRAINT_PASS;
  if (on && (fabs(d) <= tolerance)) return R3_SURFEL_CONSTRAINT_PASS;
  if (below && (d < tolerance)) return R3_SURFEL_CONSTRAINT_PASS;
  return R3_SURFEL_CONSTRAINT_FAIL;
}



////////////////////////////////////////////////////////////////////////
// GRID CONSTRAINT FUNCTIONS
////////////////////////////////////////////////////////////////////////

R3SurfelGridConstraint::
R3SurfelGridConstraint(const R3Grid *grid, 
  int comparison_type, 
  int surfel_value_type, int grid_value_type, 
  RNScalar surfel_operand, RNScalar grid_operand,
  RNScalar epsilon)
  : grid(grid),
    comparison_type(comparison_type),
    surfel_value_type(surfel_value_type),
    grid_value_type(grid_value_type),
    surfel_operand(surfel_operand),
    grid_operand(grid_operand),
    epsilon(epsilon)
{
}



int R3SurfelGridConstraint::
Check(const R3Box& box) const
{
  // Return whether any point in box can satisfy constraint
  if (R3Intersects(box, grid->WorldBox())) return R3_SURFEL_CONSTRAINT_MAYBE;
  return R3_SURFEL_CONSTRAINT_FAIL;
}



int R3SurfelGridConstraint::
Check(const R3Point& point) const
{
  // Get grid value
  RNScalar grid_value = grid_operand;
  if (grid_value_type != R3_SURFEL_CONSTRAINT_OPERAND) grid_value = grid->WorldValue(point);
  if (grid_value == R2_GRID_UNKNOWN_VALUE) return R3_SURFEL_CONSTRAINT_FAIL;

  // Get surfel value
  RNScalar surfel_value = surfel_operand;
  if (surfel_value_type == R3_SURFEL_CONSTRAINT_X) surfel_value = point.X();
  else if (surfel_value_type == R3_SURFEL_CONSTRAINT_Y) surfel_value = point.Y();
  else if (surfel_value_type == R3_SURFEL_CONSTRAINT_Z) surfel_value = point.Z();

  // Compare values
  switch (comparison_type) {
  case R3_SURFEL_CONSTRAINT_NOT_EQUAL: 
    return (RNIsNotEqual(surfel_value, grid_value, epsilon)) ?
      R3_SURFEL_CONSTRAINT_PASS : R3_SURFEL_CONSTRAINT_FAIL;
  case R3_SURFEL_CONSTRAINT_EQUAL: 
    return (RNIsEqual(surfel_value, grid_value, epsilon)) ?
      R3_SURFEL_CONSTRAINT_PASS : R3_SURFEL_CONSTRAINT_FAIL;
  case R3_SURFEL_CONSTRAINT_GREATER: 
    return (RNIsGreater(surfel_value, grid_value, epsilon)) ?
      R3_SURFEL_CONSTRAINT_PASS : R3_SURFEL_CONSTRAINT_FAIL;
  case R3_SURFEL_CONSTRAINT_GREATER_OR_EQUAL: 
    return (RNIsGreaterOrEqual(surfel_value, grid_value, epsilon)) ?
      R3_SURFEL_CONSTRAINT_PASS : R3_SURFEL_CONSTRAINT_FAIL;
  case R3_SURFEL_CONSTRAINT_LESS: 
    return (RNIsLess(surfel_value, grid_value, epsilon)) ?
      R3_SURFEL_CONSTRAINT_PASS : R3_SURFEL_CONSTRAINT_FAIL;
  case R3_SURFEL_CONSTRAINT_LESS_OR_EQUAL: 
    return (RNIsLessOrEqual(surfel_value, grid_value, epsilon)) ?
      R3_SURFEL_CONSTRAINT_PASS : R3_SURFEL_CONSTRAINT_FAIL;
  }

  // Did not pass test
  return R3_SURFEL_CONSTRAINT_FAIL;
}



////////////////////////////////////////////////////////////////////////
// PLANAR GRID CONSTRAINT FUNCTIONS
////////////////////////////////////////////////////////////////////////

R3SurfelPlanarGridConstraint::
R3SurfelPlanarGridConstraint(const R3PlanarGrid *grid, 
  RNLength max_offplane_distance, 
  int comparison_type, 
  int surfel_value_type, int grid_value_type, 
  RNScalar surfel_operand, RNScalar grid_operand,
  RNScalar epsilon)
  : grid(grid),
    max_offplane_distance(max_offplane_distance),
    comparison_type(comparison_type),
    surfel_value_type(surfel_value_type),
    grid_value_type(grid_value_type),
    surfel_operand(surfel_operand),
    grid_operand(grid_operand),
    epsilon(epsilon)
{
}



int R3SurfelPlanarGridConstraint::
Check(const R3Box& box) const
{
  // Check distance from box to plane
  RNLength offplane_distance = R3Distance(grid->Plane(), box);
  if (offplane_distance > max_offplane_distance) return R3_SURFEL_CONSTRAINT_FAIL;

  // Check distance from box to grid
  RNLength bbox_distance = R3Distance(box, grid->WorldBox());
  if (bbox_distance > max_offplane_distance) return R3_SURFEL_CONSTRAINT_FAIL;

  // Points might pass
  return R3_SURFEL_CONSTRAINT_MAYBE;
}



int R3SurfelPlanarGridConstraint::
Check(const R3Point& point) const
{
  // Check distance from point to plane
  RNLength offplane_distance = R3Distance(grid->Plane(), point);
  if (offplane_distance > max_offplane_distance) return R3_SURFEL_CONSTRAINT_FAIL;

  // Get grid value
  RNScalar grid_value = grid_operand;
  if (grid_value_type != R3_SURFEL_CONSTRAINT_OPERAND) grid_value = grid->WorldValue(point);
  if (grid_value == R2_GRID_UNKNOWN_VALUE) return R3_SURFEL_CONSTRAINT_FAIL;

  // Get surfel value
  RNScalar surfel_value = surfel_operand;
  if (surfel_value_type == R3_SURFEL_CONSTRAINT_X) surfel_value = point.X();
  else if (surfel_value_type == R3_SURFEL_CONSTRAINT_Y) surfel_value = point.Y();
  else if (surfel_value_type == R3_SURFEL_CONSTRAINT_Z) surfel_value = point.Z();

  // Compare values
  switch (comparison_type) {
  case R3_SURFEL_CONSTRAINT_NOT_EQUAL: 
    return (RNIsNotEqual(surfel_value, grid_value, epsilon)) ?
      R3_SURFEL_CONSTRAINT_PASS : R3_SURFEL_CONSTRAINT_FAIL;
  case R3_SURFEL_CONSTRAINT_EQUAL: 
    return (RNIsEqual(surfel_value, grid_value, epsilon)) ?
      R3_SURFEL_CONSTRAINT_PASS : R3_SURFEL_CONSTRAINT_FAIL;
  case R3_SURFEL_CONSTRAINT_GREATER: 
    return (RNIsGreater(surfel_value, grid_value, epsilon)) ?
      R3_SURFEL_CONSTRAINT_PASS : R3_SURFEL_CONSTRAINT_FAIL;
  case R3_SURFEL_CONSTRAINT_GREATER_OR_EQUAL: 
    return (RNIsGreaterOrEqual(surfel_value, grid_value, epsilon)) ?
      R3_SURFEL_CONSTRAINT_PASS : R3_SURFEL_CONSTRAINT_FAIL;
  case R3_SURFEL_CONSTRAINT_LESS: 
    return (RNIsLess(surfel_value, grid_value, epsilon)) ?
      R3_SURFEL_CONSTRAINT_PASS : R3_SURFEL_CONSTRAINT_FAIL;
  case R3_SURFEL_CONSTRAINT_LESS_OR_EQUAL: 
    return (RNIsLessOrEqual(surfel_value, grid_value, epsilon)) ?
      R3_SURFEL_CONSTRAINT_PASS : R3_SURFEL_CONSTRAINT_FAIL;
  }

  // Did not pass test
  return R3_SURFEL_CONSTRAINT_FAIL;
}



////////////////////////////////////////////////////////////////////////
// OVERHEAD GRID CONSTRAINT FUNCTIONS
////////////////////////////////////////////////////////////////////////

R3SurfelOverheadGridConstraint::
R3SurfelOverheadGridConstraint(const R2Grid *grid, 
  int comparison_type, 
  int surfel_value_type, int grid_value_type, 
  RNScalar surfel_operand, RNScalar grid_operand,
  RNScalar epsilon)
  : grid(grid),
    comparison_type(comparison_type),
    surfel_value_type(surfel_value_type),
    grid_value_type(grid_value_type),
    surfel_operand(surfel_operand),
    grid_operand(grid_operand),
    epsilon(epsilon)
{
}



int R3SurfelOverheadGridConstraint::
Check(const R3Box& box) const
{
  // Return whether any point in box can satisfy constraint
  R2Box projected_box(box[0][0], box[0][1], box[1][0], box[1][1]);
  if (R2Intersects(projected_box, grid->WorldBox())) return R3_SURFEL_CONSTRAINT_MAYBE;
  return R3_SURFEL_CONSTRAINT_FAIL;
}



int R3SurfelOverheadGridConstraint::
Check(const R3Point& point) const
{
  // Get grid value
  RNScalar grid_value = grid_operand;
  if (grid_value_type != R3_SURFEL_CONSTRAINT_OPERAND) grid_value = grid->WorldValue(point[0], point[1]);
  if (grid_value == R2_GRID_UNKNOWN_VALUE) return R3_SURFEL_CONSTRAINT_FAIL;

  // Get surfel value
  RNScalar surfel_value = surfel_operand;
  if (surfel_value_type == R3_SURFEL_CONSTRAINT_X) surfel_value = point.X();
  else if (surfel_value_type == R3_SURFEL_CONSTRAINT_Y) surfel_value = point.Y();
  else if (surfel_value_type == R3_SURFEL_CONSTRAINT_Z) surfel_value = point.Z();

  // Compare values
  switch (comparison_type) {
  case R3_SURFEL_CONSTRAINT_NOT_EQUAL: 
    return (RNIsNotEqual(surfel_value, grid_value, epsilon)) ?
      R3_SURFEL_CONSTRAINT_PASS : R3_SURFEL_CONSTRAINT_FAIL;
  case R3_SURFEL_CONSTRAINT_EQUAL: 
    return (RNIsEqual(surfel_value, grid_value, epsilon)) ?
      R3_SURFEL_CONSTRAINT_PASS : R3_SURFEL_CONSTRAINT_FAIL;
  case R3_SURFEL_CONSTRAINT_GREATER: 
    return (RNIsGreater(surfel_value, grid_value, epsilon)) ?
      R3_SURFEL_CONSTRAINT_PASS : R3_SURFEL_CONSTRAINT_FAIL;
  case R3_SURFEL_CONSTRAINT_GREATER_OR_EQUAL: 
    return (RNIsGreaterOrEqual(surfel_value, grid_value, epsilon)) ?
      R3_SURFEL_CONSTRAINT_PASS : R3_SURFEL_CONSTRAINT_FAIL;
  case R3_SURFEL_CONSTRAINT_LESS: 
    return (RNIsLess(surfel_value, grid_value, epsilon)) ?
      R3_SURFEL_CONSTRAINT_PASS : R3_SURFEL_CONSTRAINT_FAIL;
  case R3_SURFEL_CONSTRAINT_LESS_OR_EQUAL: 
    return (RNIsLessOrEqual(surfel_value, grid_value, epsilon)) ?
      R3_SURFEL_CONSTRAINT_PASS : R3_SURFEL_CONSTRAINT_FAIL;
  }

  // Did not pass test
  return R3_SURFEL_CONSTRAINT_FAIL;
}



////////////////////////////////////////////////////////////////////////
// VIEW CONSTRAINT FUNCTIONS
////////////////////////////////////////////////////////////////////////

R3SurfelViewConstraint::
R3SurfelViewConstraint(
  const R3Viewer& viewer,
  const R2Grid *image_mask,
  RNBoolean must_be_inside_frustum)
  : world_to_camera(viewer.Camera().CoordSystem().InverseMatrix()),
    w(viewer.Viewport().Width()),
    h(viewer.Viewport().Height()),
    fx(0),
    fy(0),
    cx(0.5*viewer.Viewport().Width()),
    cy(0.5*viewer.Viewport().Height()),
    neardist(viewer.Camera().Near()),
    fardist(viewer.Camera().Far()),
    image_mask(image_mask),
    must_be_inside_frustum(must_be_inside_frustum)
{
  // Get convenient variables
  const R3Camera& camera = viewer.Camera();
  R3Point world_viewpoint = camera.Origin();
  R3Vector world_towards = camera.Towards();
  R3Vector world_up = camera.Up();
  R3Vector world_right = camera.Right();
  RNScalar xfov = camera.XFOV();
  RNScalar yfov = camera.YFOV();
  RNScalar tan_xfov = tan(xfov);
  RNScalar tan_yfov = tan(yfov);

  // Update focal lengths
  this->fx = (tan_xfov != 0) ? 0.5 * w / tan(xfov) : 0;
  this->fy = (tan_yfov != 0) ? 0.5 * h / tan(yfov) : 0;

  // Compute frustum halfspaces
  R3Vector world_normal;
  world_normal = world_right;  world_normal.Rotate(world_up, xfov);
  frustum[RN_LO][RN_X] = R3Halfspace(world_viewpoint, world_normal);
  world_normal = -world_right;  world_normal.Rotate(world_up, -xfov);
  frustum[RN_HI][RN_X] = R3Halfspace(world_viewpoint, world_normal);
  world_normal = world_up;  world_normal.Rotate(world_right, -yfov);
  frustum[RN_LO][RN_Y] = R3Halfspace(world_viewpoint, world_normal);
  world_normal = -world_up;  world_normal.Rotate(world_right, yfov);
  frustum[RN_HI][RN_Y] = R3Halfspace(world_viewpoint, world_normal);
  frustum[RN_LO][RN_Z] = R3Halfspace(world_viewpoint + neardist*world_towards, world_towards);
  frustum[RN_HI][RN_Z] = R3Halfspace(world_viewpoint + fardist*world_towards, -world_towards);
}



R3SurfelViewConstraint::
R3SurfelViewConstraint(
  const R3SurfelImage& image,
  RNLength neardist, RNLength fardist, 
  const R2Grid *image_mask,
  RNBoolean must_be_inside_frustum)
  : world_to_camera(image.Extrinsics()),
    w(image.ImageWidth()), h(image.ImageHeight()),
    fx(image.XFocal()), fy(image.YFocal()),
    cx(image.ImageCenter().X()), cy(image.ImageCenter().Y()),
    neardist(neardist), fardist(fardist),
    image_mask(image_mask),
    must_be_inside_frustum(must_be_inside_frustum)
{
  // Get convenient variables
  R3Point world_viewpoint = image.Viewpoint();
  R3Vector world_towards = image.Towards();
  R3Vector world_up = image.Up();
  R3Vector world_right = image.Right();
  RNScalar xfov = image.XFOV();
  RNScalar yfov = image.YFOV();

  // Compute default image center
  if (this->cx <= 0) this->cx = 0.5*w;
  if (this->cy <= 0) this->cy = 0.5*h;

  // Compute default near and far distances
  if (this->neardist <= 0) this->neardist = RN_EPSILON;
  if (this->fardist <= 0) this->fardist = RN_INFINITY;

  // Compute frustum halfspaces
  R3Vector world_normal;
  world_normal = world_right;  world_normal.Rotate(world_up, xfov);
  frustum[RN_LO][RN_X] = R3Halfspace(world_viewpoint, world_normal);
  world_normal = -world_right;  world_normal.Rotate(world_up, -xfov);
  frustum[RN_HI][RN_X] = R3Halfspace(world_viewpoint, world_normal);
  world_normal = world_up;  world_normal.Rotate(world_right, -yfov);
  frustum[RN_LO][RN_Y] = R3Halfspace(world_viewpoint, world_normal);
  world_normal = -world_up;  world_normal.Rotate(world_right, yfov);
  frustum[RN_HI][RN_Y] = R3Halfspace(world_viewpoint, world_normal);
  frustum[RN_LO][RN_Z] = R3Halfspace(world_viewpoint + neardist*world_towards, world_towards);
  frustum[RN_HI][RN_Z] = R3Halfspace(world_viewpoint + fardist*world_towards, -world_towards);
}



R3SurfelViewConstraint::
R3SurfelViewConstraint(const R3Point& world_viewpoint,
  const R3Vector& world_towards, const R3Vector& world_up,
  int w, int h, RNScalar fx, RNScalar fy, RNScalar cx, RNScalar cy,
  RNLength neardist, RNLength fardist, 
  const R2Grid *image_mask,
  RNBoolean must_be_inside_frustum)
  : world_to_camera(),
    w(w), h(h),
    fx(fx), fy(fy),
    cx(cx), cy(cy),
    neardist(neardist), fardist(fardist),
    image_mask(image_mask),
    must_be_inside_frustum(must_be_inside_frustum)
{
  // Compute default image center
  if (this->cx <= 0) this->cx = 0.5*w;
  if (this->cy <= 0) this->cy = 0.5*h;

  // Compute default near and far distances
  if (this->neardist <= 0) this->neardist = RN_EPSILON;
  if (this->fardist <= 0) this->fardist = RN_INFINITY;

  // Compute world_to_camera matrix
  R3CoordSystem cs(world_viewpoint, R3Triad(world_towards, world_up));
  world_to_camera = cs.InverseMatrix();
  
  // Compute frustum halfspaces
  R3Vector world_normal;
  RNScalar xfov = atan(0.5*w/fx);
  RNScalar yfov = atan(0.5*h/fy);
  R3Vector world_right = world_towards % world_up; world_right.Normalize();
  world_normal = world_right;  world_normal.Rotate(world_up, xfov);
  frustum[RN_LO][RN_X] = R3Halfspace(world_viewpoint, world_normal);
  world_normal = -world_right;  world_normal.Rotate(world_up, -xfov);
  frustum[RN_HI][RN_X] = R3Halfspace(world_viewpoint, world_normal);
  world_normal = world_up;  world_normal.Rotate(world_right, -yfov);
  frustum[RN_LO][RN_Y] = R3Halfspace(world_viewpoint, world_normal);
  world_normal = -world_up;  world_normal.Rotate(world_right, yfov);
  frustum[RN_HI][RN_Y] = R3Halfspace(world_viewpoint, world_normal);
  frustum[RN_LO][RN_Z] = R3Halfspace(world_viewpoint + neardist*world_towards, world_towards);
  frustum[RN_HI][RN_Z] = R3Halfspace(world_viewpoint + fardist*world_towards, -world_towards);
}



int R3SurfelViewConstraint::
Check(const R3Box& box) const
{
  // Return whether any point in box can satisfy constraint (not implemented yet)
  if (must_be_inside_frustum) {
    for (int dir = 0; dir < 2; dir++) {
      for (int dim = 0; dim < 3; dim++) {
        if (!R3Intersects(frustum[dir][dim], box)) {
          return R3_SURFEL_CONSTRAINT_FAIL;
        }
      }
    }
  }
  
  // Passed all tests
  return R3_SURFEL_CONSTRAINT_MAYBE;
}



int R3SurfelViewConstraint::
Check(const R3Point& point) const
{
  // Transform into camera coordinates (z is negative)
  R3Point camera_point = world_to_camera * point;
  RNScalar depth = -camera_point.Z();

  // Check if z is too near or far
  if (depth < neardist) return R3_SURFEL_CONSTRAINT_FAIL;
  if (depth > fardist) return R3_SURFEL_CONSTRAINT_FAIL;

  // Compute image position
  int ix = (int) (cx + fx*camera_point.X()/depth);
  int iy = (int) (cy + fy*camera_point.Y()/depth);

  // Check image position
  if (must_be_inside_frustum) {
    // Return fail if outside
    if ((ix < 0) || (ix >= w)) return R3_SURFEL_CONSTRAINT_FAIL;
    if ((iy < 0) || (iy >= h)) return R3_SURFEL_CONSTRAINT_FAIL;
  }
  else {
    // Snap to closest pixel on boundary
    if (ix < 0) ix = 0;
    if (ix >= w) ix = w-1;
    if (iy < 0) iy = 0;
    if (iy >= h) iy = h-1;
  }
  
  // Check if point is masked
  if (image_mask) {
    RNScalar mask_value = image_mask->GridValue(ix, iy);
    if (mask_value == R2_GRID_UNKNOWN_VALUE) return R3_SURFEL_CONSTRAINT_FAIL;
    if (mask_value == 0.0) return R3_SURFEL_CONSTRAINT_FAIL;
  }
  
  // Passed all tests
  return R3_SURFEL_CONSTRAINT_PASS;
}



////////////////////////////////////////////////////////////////////////
// MESH CONSTRAINT FUNCTIONS
////////////////////////////////////////////////////////////////////////

R3SurfelMeshConstraint::
R3SurfelMeshConstraint(R3Mesh *mesh, const R3Affine& surfels_to_mesh, RNLength max_distance)
  : mesh(mesh),
    surfels_to_mesh(surfels_to_mesh),
    max_distance(max_distance),
    tree(NULL)
{
  // Just checking
  assert(mesh);
  tree = new R3MeshSearchTree(mesh);
  assert(tree);
}


R3SurfelMeshConstraint::
~R3SurfelMeshConstraint(void)
{
  // Delete search tree
  if (tree) delete tree;
}



int R3SurfelMeshConstraint::
Check(const R3Box& box) const
{
  // Return whether any point in box can satisfy constraint
  R3Box transformed_box = box;
  transformed_box.Transform(surfels_to_mesh);
  RNLength distance = R3Distance(transformed_box, mesh->BBox());
  if (distance > max_distance) return R3_SURFEL_CONSTRAINT_FAIL;
  return R3_SURFEL_CONSTRAINT_MAYBE;
}



int R3SurfelMeshConstraint::
Check(const R3Point& point) const
{
  // Check mesh
  if (!mesh) return R3_SURFEL_CONSTRAINT_FAIL;

  // Transform position
  R3Point transformed_point = point;
  transformed_point.Transform(surfels_to_mesh);

  // Check if closest point is within max_distance
  R3MeshIntersection closest;
  tree->FindClosest(transformed_point, closest, 0, max_distance);
  if (closest.type == R3_MESH_NULL_TYPE) return R3_SURFEL_CONSTRAINT_FAIL;
  return R3_SURFEL_CONSTRAINT_PASS;
}



////////////////////////////////////////////////////////////////////////
// OBJECT CONSTRAINT FUNCTIONS
////////////////////////////////////////////////////////////////////////

R3SurfelObjectConstraint::
R3SurfelObjectConstraint(R3SurfelObject *target_object, RNBoolean converse)
  : target_object(target_object),
    converse(converse)
{
}



int R3SurfelObjectConstraint::
Check(const R3SurfelObject *object) const
{
  // Return whether surfel satisfies constraint
  int status = R3_SURFEL_CONSTRAINT_PASS;
  if (!object) status = R3_SURFEL_CONSTRAINT_FAIL;
  else if (target_object && (object != target_object)) status = R3_SURFEL_CONSTRAINT_FAIL;
  if (converse) {
    if (status == R3_SURFEL_CONSTRAINT_FAIL) status = R3_SURFEL_CONSTRAINT_PASS;
    else status = R3_SURFEL_CONSTRAINT_FAIL;
  }
  return status;
}



int R3SurfelObjectConstraint::
Check(const R3SurfelNode *node) const
{
  // Return whether surfel satisfies constraint
  return Check(node->Object());
}



int R3SurfelObjectConstraint::
Check(const R3SurfelBlock *block) const
{
  // Return whether surfel satisfies constraint
  R3SurfelNode *node = block->Node();
  R3SurfelObject *object = (node) ? node->Object() : NULL;
  return Check(object);
}



int R3SurfelObjectConstraint::
Check(const R3SurfelBlock *block, const R3Surfel *surfel) const
{
  // Return whether surfel satisfies constraint
  return R3_SURFEL_CONSTRAINT_PASS;
}



////////////////////////////////////////////////////////////////////////
// LABEL CONSTRAINT FUNCTIONS
////////////////////////////////////////////////////////////////////////

R3SurfelLabelConstraint::
R3SurfelLabelConstraint(R3SurfelLabel *target_label, RNBoolean converse)
  : target_label(target_label),
    converse(converse)
{
}



int R3SurfelLabelConstraint::
Check(const R3SurfelObject *object) const
{
  // Return whether surfel satisfies constraint
  int status = R3_SURFEL_CONSTRAINT_PASS;
  if (!object) status = R3_SURFEL_CONSTRAINT_FAIL;
  else {
    R3SurfelLabel *label = object->CurrentLabel();
    if (!label) status = R3_SURFEL_CONSTRAINT_FAIL;
    else {
      if (target_label && (label != target_label)) status = R3_SURFEL_CONSTRAINT_FAIL;
    }
  }
  if (converse) {
    if (status == R3_SURFEL_CONSTRAINT_FAIL) status = R3_SURFEL_CONSTRAINT_PASS;
    else status = R3_SURFEL_CONSTRAINT_FAIL;
  }
  return status;
}



int R3SurfelLabelConstraint::
Check(const R3SurfelNode *node) const
{
  // Return whether surfel satisfies constraint
  return Check(node->Object());
}



int R3SurfelLabelConstraint::
Check(const R3SurfelBlock *block) const
{
  // Return whether surfel satisfies constraint
  R3SurfelNode *node = block->Node();
  R3SurfelObject *object = (node) ? node->Object() : NULL;
  return Check(object);
}



int R3SurfelLabelConstraint::
Check(const R3SurfelBlock *block, const R3Surfel *surfel) const
{
  // Return whether surfel satisfies constraint
  return R3_SURFEL_CONSTRAINT_PASS;
}



////////////////////////////////////////////////////////////////////////
// BOUNDARY CONSTRAINT FUNCTIONS
////////////////////////////////////////////////////////////////////////

R3SurfelBoundaryConstraint::
R3SurfelBoundaryConstraint(RNBoolean include_interior, RNBoolean include_border,
  RNBoolean include_silhouette, RNBoolean include_shadow)
  : include_interior(include_interior), include_border(include_border),
    include_silhouette(include_silhouette), include_shadow(include_shadow)
{
}



int R3SurfelBoundaryConstraint::
Check(const R3SurfelBlock *block, const R3Surfel *surfel) const
{
  // Return whether surfel satisfies constraint
  if (include_interior && !(surfel->Flags() & R3_SURFEL_BOUNDARY_FLAGS)) return R3_SURFEL_CONSTRAINT_PASS;
  if (include_silhouette && surfel->IsOnSilhouetteBoundary()) return R3_SURFEL_CONSTRAINT_PASS;
  if (include_shadow && surfel->IsOnShadowBoundary()) return R3_SURFEL_CONSTRAINT_PASS;
  if (include_border && surfel->IsOnBorderBoundary()) return R3_SURFEL_CONSTRAINT_PASS;
  return R3_SURFEL_CONSTRAINT_FAIL;
}



////////////////////////////////////////////////////////////////////////
// SOURCE CONSTRAINT FUNCTIONS
////////////////////////////////////////////////////////////////////////

R3SurfelSourceConstraint::
R3SurfelSourceConstraint(RNBoolean include_aerial, RNBoolean include_terrestrial)
  : include_aerial(include_aerial), include_terrestrial(include_terrestrial)
{
}



int R3SurfelSourceConstraint::
Check(const R3SurfelNode *node) const
{
  // Return whether surfel satisfies constraint
  if (include_terrestrial && node->HasTerrestrial()) return R3_SURFEL_CONSTRAINT_MAYBE;
  if (include_aerial && node->HasAerial()) return R3_SURFEL_CONSTRAINT_MAYBE;
  return R3_SURFEL_CONSTRAINT_FAIL;
}



int R3SurfelSourceConstraint::
Check(const R3SurfelBlock *block) const
{
  // Return whether surfel satisfies constraint
  if (include_terrestrial && block->HasTerrestrial()) return R3_SURFEL_CONSTRAINT_MAYBE;
  if (include_aerial && block->HasAerial()) return R3_SURFEL_CONSTRAINT_MAYBE;
  return R3_SURFEL_CONSTRAINT_FAIL;
}



int R3SurfelSourceConstraint::
Check(const R3SurfelBlock *block, const R3Surfel *surfel) const
{
  // Return whether surfel satisfies constraint
  if (include_terrestrial && surfel->IsTerrestrial()) return R3_SURFEL_CONSTRAINT_PASS;
  if (include_aerial && surfel->IsAerial()) return R3_SURFEL_CONSTRAINT_PASS;
  return R3_SURFEL_CONSTRAINT_FAIL;
}



////////////////////////////////////////////////////////////////////////
// MARK CONSTRAINT FUNCTIONS
////////////////////////////////////////////////////////////////////////

R3SurfelMarkConstraint::
R3SurfelMarkConstraint(RNBoolean include_marked, RNBoolean include_unmarked)
  : include_marked(include_marked), include_unmarked(include_unmarked)
{
}



int R3SurfelMarkConstraint::
Check(const R3SurfelBlock *block) const
{
  // Get block database
  R3SurfelDatabase *database = block->Database();
  if (!database) return R3_SURFEL_CONSTRAINT_MAYBE;

  // Check whether block is resident
  RNBoolean resident = database->IsBlockResident((R3SurfelBlock *) block);
  if (!resident) {
    if (include_unmarked) return R3_SURFEL_CONSTRAINT_PASS; 
    else return R3_SURFEL_CONSTRAINT_FAIL; 
  }

  // Block is resident -- check surfels
  return R3_SURFEL_CONSTRAINT_MAYBE;
}



int R3SurfelMarkConstraint::
Check(const R3SurfelBlock *block, const R3Surfel *surfel) const
{
  // Return whether surfel satisfies constraint
  if (surfel->IsMarked()) {
    if (include_marked) return R3_SURFEL_CONSTRAINT_PASS; 
  }
  else {
    if (include_unmarked) return R3_SURFEL_CONSTRAINT_PASS;
  }
  return R3_SURFEL_CONSTRAINT_FAIL;
}



////////////////////////////////////////////////////////////////////////
// MULTI CONSTRAINT FUNCTIONS
////////////////////////////////////////////////////////////////////////

R3SurfelMultiConstraint::
R3SurfelMultiConstraint(void)
  : constraints()
{
}



void R3SurfelMultiConstraint::
InsertConstraint(const R3SurfelConstraint *constraint) 
{
  // Insert constraint
  constraints.Insert(constraint);
}



void R3SurfelMultiConstraint::
RemoveConstraint(const R3SurfelConstraint *constraint) 
{
  // Remove constraint
  constraints.Remove(constraint);
}




int R3SurfelMultiConstraint::
Check(const R3SurfelObject *object) const
{
  // Initialize return status 
  int status = R3_SURFEL_CONSTRAINT_PASS;

  // Check whether all constraints are met
  for (int i = 0; i < constraints.NEntries(); i++) {
    const R3SurfelConstraint *constraint = constraints.Kth(i);
    int s = constraint->Check(object);
    if (s == R3_SURFEL_CONSTRAINT_FAIL) return R3_SURFEL_CONSTRAINT_FAIL;
    else if (s == R3_SURFEL_CONSTRAINT_MAYBE) status = R3_SURFEL_CONSTRAINT_MAYBE;
  }

  // Return how all constraints are met
  return status;
}



int R3SurfelMultiConstraint::
Check(const R3SurfelNode *node) const
{
  // Initialize return status 
  int status = R3_SURFEL_CONSTRAINT_PASS;

  // Check whether all constraints are met
  for (int i = 0; i < constraints.NEntries(); i++) {
    const R3SurfelConstraint *constraint = constraints.Kth(i);
    int s = constraint->Check(node);
    if (s == R3_SURFEL_CONSTRAINT_FAIL) return R3_SURFEL_CONSTRAINT_FAIL;
    else if (s == R3_SURFEL_CONSTRAINT_MAYBE) status = R3_SURFEL_CONSTRAINT_MAYBE;
  }

  // Return how all constraints are met
  return status;
}



int R3SurfelMultiConstraint::
Check(const R3SurfelBlock *block) const
{
  // Initialize return status 
  int status = R3_SURFEL_CONSTRAINT_PASS;

  // Check whether all constraints are met
  for (int i = 0; i < constraints.NEntries(); i++) {
    const R3SurfelConstraint *constraint = constraints.Kth(i);
    int s = constraint->Check(block);
    if (s == R3_SURFEL_CONSTRAINT_FAIL) return R3_SURFEL_CONSTRAINT_FAIL;
    else if (s == R3_SURFEL_CONSTRAINT_MAYBE) status = R3_SURFEL_CONSTRAINT_MAYBE;
  }

  // Return how all constraints are met
  return status;
}



int R3SurfelMultiConstraint::
Check(const R3SurfelBlock *block, const R3Surfel *surfel) const
{
  // Check whether all constraints are met
  for (int i = 0; i < constraints.NEntries(); i++) {
    const R3SurfelConstraint *constraint = constraints.Kth(i);
    if (!constraint->Check(block, surfel)) return R3_SURFEL_CONSTRAINT_FAIL;
  }

  // Return how all constraints are met
  return R3_SURFEL_CONSTRAINT_PASS;
}



int R3SurfelMultiConstraint::
Check(const R3Box& box) const
{
  // Initialize return status 
  int status = R3_SURFEL_CONSTRAINT_PASS;

  // Check whether all constraints are met
  for (int i = 0; i < constraints.NEntries(); i++) {
    const R3SurfelConstraint *constraint = constraints.Kth(i);
    int s = constraint->Check(box);
    if (s == R3_SURFEL_CONSTRAINT_FAIL) return R3_SURFEL_CONSTRAINT_FAIL;
    else if (s == R3_SURFEL_CONSTRAINT_MAYBE) status = R3_SURFEL_CONSTRAINT_MAYBE;
  }

  // Return how all constraints are met
  return status;
}



int R3SurfelMultiConstraint::
Check(const R3Point& point) const
{
  // Check whether all constraints are met
  for (int i = 0; i < constraints.NEntries(); i++) {
    const R3SurfelConstraint *constraint = constraints.Kth(i);
    if (!constraint->Check(point)) return R3_SURFEL_CONSTRAINT_FAIL;
  }

  // All constraints are met
  return R3_SURFEL_CONSTRAINT_PASS;
}



} // namespace gaps
