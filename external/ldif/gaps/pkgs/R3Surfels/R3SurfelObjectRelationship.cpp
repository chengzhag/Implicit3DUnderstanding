/* Source file for the R3 surfel object relationship class */



////////////////////////////////////////////////////////////////////////
// INCLUDE FILES
////////////////////////////////////////////////////////////////////////

#include "R3Surfels.h"



////////////////////////////////////////////////////////////////////////
// Namespace
////////////////////////////////////////////////////////////////////////

namespace gaps {



////////////////////////////////////////////////////////////////////////
// CONSTRUCTORS/DESTRUCTORS
////////////////////////////////////////////////////////////////////////

R3SurfelObjectRelationship::
R3SurfelObjectRelationship(int type, const RNArray<R3SurfelObject *>& objects, RNScalar *operands, int noperands)
  : scene(NULL),
    scene_index(-1),
    objects(objects),
    operands(NULL),
    noperands(0),
    type(type)
{
  // Check if operands were provided
  if ((noperands > 0) && (operands)) {
    // Copy operands
    this->noperands = noperands;
    this->operands = new RNScalar [ this->noperands ];
    for (int i = 0; i < this->noperands; i++) {
      this->operands[i] = operands[i];
    }
  }
  else {
    // Compute operands
    UpdateOperands();
  }
}



R3SurfelObjectRelationship::
R3SurfelObjectRelationship(int type, R3SurfelObject *object0, R3SurfelObject *object1, RNScalar *operands, int noperands)
  : scene(NULL),
    scene_index(-1),
    objects(),
    operands(NULL),
    noperands(0),
    type(type)
{
  // Insert objects
  objects.Insert(object0);
  objects.Insert(object1);

  // Check if operands were provided
  if ((noperands > 0) && (operands)) {
    // Copy operands
    this->noperands = noperands;
    this->operands = new RNScalar [ this->noperands ];
    for (int i = 0; i < this->noperands; i++) {
      this->operands[i] = operands[i];
    }
  }
  else {
    // Compute operands
    UpdateOperands();
  }
}



R3SurfelObjectRelationship::
R3SurfelObjectRelationship(const R3SurfelObjectRelationship& relationship)
  : scene(NULL),
    scene_index(-1),
    objects(relationship.objects),
    operands(NULL),
    noperands(0),
    type(relationship.type)
{
  // Copy operands
  if ((relationship.noperands > 0) && (relationship.operands)) {
    this->noperands = relationship.noperands;
    this->operands = new RNScalar [ this->noperands ];
    for (int i = 0; i < this->noperands; i++) {
      this->operands[i] = relationship.operands[i];
    }
  }
}



R3SurfelObjectRelationship::
~R3SurfelObjectRelationship(void)
{
  // Remove from scene (and everything in it)
  if (scene) scene->RemoveObjectRelationship(this);

  // Delete operands
  if (operands) delete [] operands;
}






////////////////////////////////////////////////////////////////////////
// DISPLAY FUNCDTIONS
////////////////////////////////////////////////////////////////////////

void R3SurfelObjectRelationship::
Draw(RNFlags flags) const
{
  // Draw property based on type
  switch (type) {
  case R3_SURFEL_OBJECT_OVERLAP_RELATIONSHIP: 
    // Draw line between objects
    if (objects.NEntries() == 2) {
      if (flags[R3_SURFEL_COLOR_DRAW_FLAG]) glColor3d(0, 1, 0);
      glBegin(GL_LINES);
      R3LoadPoint(objects[0]->Centroid());
      R3LoadPoint(objects[1]->Centroid());
      glEnd();
      break;
    }
  }
}



////////////////////////////////////////////////////////////////////////
// INTERNAL FUNCDTIONS
////////////////////////////////////////////////////////////////////////

void R3SurfelObjectRelationship::
UpdateOperands(void) 
{
  switch (type) {
  }
}



} // namespace gaps
  
