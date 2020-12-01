/* Source file for the R3 surfel class */



/* Include files */

#include "R3Surfels.h"



// Namespace

namespace gaps {



/* Public functions */

R3Surfel::
R3Surfel(void)
  : timestamp(0),
    identifier(0),
    flags(0)
{
  // Set everything
  this->position[0] = 0;
  this->position[1] = 0;
  this->position[2] = 0;
  this->normal[0] = 0;
  this->normal[1] = 0;
  this->normal[2] = 0;
  this->tangent[0] = 0;
  this->tangent[1] = 0;
  this->tangent[2] = 0;
  this->radius[0] = 0;
  this->radius[1] = 0;
  this->color[0] = 0;
  this->color[1] = 0;
  this->color[2] = 0;
}



R3Surfel::
R3Surfel(float x, float y, float z, 
  unsigned char r, unsigned char g, unsigned char b, 
  RNBoolean aerial)
  : timestamp(0),
    identifier(0),
    flags(0)
{
  // Set everything
  this->position[0] = x;
  this->position[1] = y;
  this->position[2] = z;
  this->normal[0] = 0;
  this->normal[1] = 0;
  this->normal[2] = 0;
  this->tangent[0] = 0;
  this->tangent[1] = 0;
  this->tangent[2] = 0;
  this->radius[0] = 0;
  this->radius[1] = 0;
  this->color[0] = r;
  this->color[1] = g;
  this->color[2] = b;
  SetAerial(aerial);
}



R3Surfel::
R3Surfel(float x, float y, float z, float nx, float ny, float nz,
  float radius, unsigned char r, unsigned char g, unsigned char b, 
  unsigned char flags)
  : timestamp(0),
    identifier(0),
    flags(flags)
{
  // Set everything
  this->position[0] = x;
  this->position[1] = y;
  this->position[2] = z;
  this->normal[0] = (RNInt16) (32767.0 * nx + 0.5);
  this->normal[1] = (RNInt16) (32767.0 * ny + 0.5);
  this->normal[2] = (RNInt16) (32767.0 * nz + 0.5);
  this->tangent[0] = 0;
  this->tangent[1] = 0;
  this->tangent[2] = 0;
  this->radius[0] = (RNUInt16) (8192.0 * radius + 0.5);
  this->radius[1] = this->radius[0];
  this->color[0] = r;
  this->color[1] = g;
  this->color[2] = b;

  // Update flags
  if ((nx != 0) && (ny != 0) && (nz != 0)) {
    this->flags |= R3_SURFEL_NORMAL_FLAG;
  }
}



R3Surfel::
R3Surfel(float x, float y, float z,
  float nx, float ny, float nz,
  float tx, float ty, float tz,
  float radius0, float radius1,
  unsigned char r, unsigned char g, unsigned char b, 
  float timestamp, unsigned int identifier,
  unsigned char flags)
  : timestamp(timestamp),
    identifier(identifier),
    flags(flags)
{
  // Set everything
  this->position[0] = x;
  this->position[1] = y;
  this->position[2] = z;
  this->normal[0] = (RNInt16) (32767.0 * nx + 0.5);
  this->normal[1] = (RNInt16) (32767.0 * ny + 0.5);
  this->normal[2] = (RNInt16) (32767.0 * nz + 0.5);
  this->tangent[0] = (RNInt16) (32767.0 * tx + 0.5);
  this->tangent[1] = (RNInt16) (32767.0 * ty + 0.5);
  this->tangent[2] = (RNInt16) (32767.0 * tz + 0.5);
  this->radius[0] = (RNUInt16) (8192.0 * radius0 + 0.5);
  this->radius[1] = (RNUInt16) (8192.0 * radius1 + 0.5);
  this->color[0] = r;
  this->color[1] = g;
  this->color[2] = b;

  // Update flags
  if ((nx != 0) && (ny != 0) && (nz != 0)) {
    this->flags |= R3_SURFEL_NORMAL_FLAG;
  }
  if ((tx != 0) && (ty != 0) && (tz != 0)) {
    this->flags |= R3_SURFEL_TANGENT_FLAG;
  }
}



} // namespace gaps
