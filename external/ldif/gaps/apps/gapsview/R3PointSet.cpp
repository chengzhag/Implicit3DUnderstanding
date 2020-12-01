// Source file for 3D point set class


namespace gaps {};
using namespace gaps;
#include "R3Shapes/R3Shapes.h"
#include "R3PointSet.h"



R3PointSet::
R3PointSet(void)
  : points(),
    normals(),
    values(),
    bbox(FLT_MAX, FLT_MAX, FLT_MAX, -FLT_MAX, -FLT_MAX, -FLT_MAX)
{
}



R3PointSet::
~R3PointSet(void)
{
}



R3Point R3PointSet::
Centroid(void) const
{
  // Return centroid of point set
  R3Point centroid = R3zero_point;
  for (unsigned int i = 0; i < points.size(); i++) 
    centroid += points[i];
  if (points.size() > 0) centroid /= points.size();
  return centroid;
}



int R3PointSet::
FindPointIndex(const R3Point& point, double tolerance) const
{
  // Return first point within tolerance
  double squared_tolerance = tolerance * tolerance;
  for (unsigned int i = 0; i < points.size(); i++) {
    if (R3SquaredDistance(points[i], point) <= squared_tolerance) {
      return i;
    }
  }

  // No point found
  return -1;
}



void R3PointSet::
RemovePoint(int index)
{
  // Remove point at given index
  points.erase(points.begin()+index);
  normals.erase(normals.begin()+index);
  values.erase(values.begin()+index);
}



void R3PointSet::
InsertPoint(const R3Point& position, const R3Vector& normal, RNScalar value)
{
  // Insert point
  points.push_back(position);
  normals.push_back(normal);
  values.push_back(value);
  bbox.Union(position);
}



void R3PointSet::
SetPointPosition(int index, const R3Point& position)
{
  // Set position for point with given index
  points[index] = position;
  bbox.Union(position);
}



void R3PointSet::
SetPointNormal(int index, const R3Vector& normal)
{
  // Set normal for point with given index
  normals[index] = normal;
}



void R3PointSet::
SetPointValue(int index, RNScalar value)
{
  // Set value for point with given index
  values[index] = value;
}



int R3PointSet::
ReadFile(const char *filename)
{
  // Parse input filename extension
  const char *extension;
  if (!(extension = strrchr(filename, '.'))) {
    printf("Filename %s has no extension (e.g., .mm)\n", filename);
    return 0;
  }

  // Read file of appropriate type
  if (!strcmp(extension, ".xyzn")) return ReadASCIIFile(filename);
  else if (!strcmp(extension, ".pts")) return ReadBinaryFile(filename);
  else if (!strcmp(extension, ".sdf")) return ReadSDFFile(filename);
  else return ReadMeshFile(filename);

  // Should not get here
  RNFail("Unrecognized file extension in %s\n", filename);
  return 0;
}


int R3PointSet::
WriteFile(const char *filename) const
{
  // Parse input filename extension
  const char *extension;
  if (!(extension = strrchr(filename, '.'))) {
    printf("Filename %s has no extension (e.g., .mm)\n", filename);
    return 0;
  }

  // Write file of appropriate type
  if (!strcmp(extension, ".xyzn")) return WriteASCIIFile(filename);
  else if (!strcmp(extension, ".pts")) return WriteBinaryFile(filename);
  else if (!strcmp(extension, ".sdf")) return WriteSDFFile(filename);
  else return WriteMeshFile(filename);

  // Should not get here
  RNFail("Unrecognized file extension in %s\n", filename);
  return 0;
}



////////////////////////////////////////////////////////////////////////

int R3PointSet::
ReadASCIIFile(const char *filename)
{
  // Open file
  FILE *fp = fopen(filename, "r");
  if (!fp) {
    RNFail("Unable to open %s\n", filename);
    return 0;
  }

  // Read points
  double px, py, pz, nx, ny, nz;
  while (fscanf(fp, "%lf%lf%lf%lf%lf%lf", &px, &py, &pz, &nx, &ny, &nz) == (unsigned int) 6) {
    R3Point position(px,py,pz);
    R3Vector normal(nx,ny,nz);
    normal.Normalize();
    InsertPoint(position, normal);
  }

  // Close file
  fclose(fp);

  // Return success
  return 1;
}



int R3PointSet::
WriteASCIIFile(const char *filename) const
{
 // Open file
  FILE *fp = fopen(filename, "w");
  if (!fp) {
    RNFail("Unable to open %s\n", filename);
    return 0;
  }

  // Write points
  for (int i = 0; i < NPoints(); i++) {
    const R3Point& p = PointPosition(i);
    const R3Vector& n = PointNormal(i);
    fprintf(fp, "%g %g %g %g %g %g\n", p[0], p[1], p[2], n[0], n[1], n[2]);
  }

  // Close file
  fclose(fp);

  // Return success
  return 1;
}



////////////////////////////////////////////////////////////////////////

int R3PointSet::
ReadBinaryFile(const char *filename)
{
  // Open file
  FILE *fp = fopen(filename, "rb");
  if (!fp) {
    RNFail("Unable to open %s\n", filename);
    return 0;
  }

  // Read points
  float buf[6];
  while (fread(buf, sizeof(float), 6, fp) == (unsigned int) 6) {
    R3Point position(buf[0], buf[1], buf[2]);
    R3Vector normal(buf[3], buf[4], buf[5]);
    normal.Normalize();
    InsertPoint(position, normal);
  }

  // Close file
  fclose(fp);

  // Return success
  return 1;
}



int R3PointSet::
WriteBinaryFile(const char *filename) const
{
 // Open file
  FILE *fp = fopen(filename, "wb");
  if (!fp) {
    RNFail("Unable to open %s\n", filename);
    return 0;
  }

  // Write points
  float buf[6];
  for (int i = 0; i < NPoints(); i++) {
    const R3Point& p = PointPosition(i);
    const R3Vector& n = PointNormal(i);
    buf[0] = p[0]; buf[1] = p[1]; buf[2] = p[2];
    buf[3] = n[0]; buf[4] = n[1]; buf[5] = n[2];
    if (fwrite(buf, sizeof(float), 6, fp) != (unsigned int) 6) {
      RNFail("Unable to write point %d to %s\n", i, filename);
      return 0;
    }
  }

  // Close file
  fclose(fp);

  // Return success
  return 1;
}



////////////////////////////////////////////////////////////////////////

int R3PointSet::
ReadSDFFile(const char *filename)
{
  // Open file
  FILE *fp = fopen(filename, "rb");
  if (!fp) {
    RNFail("Unable to open %s\n", filename);
    return 0;
  }

  // Read points
  float buf[4];
  while (fread(buf, sizeof(float), 4, fp) == (unsigned int) 4) {
    R3Point position(buf[0], buf[1], buf[2]);
    InsertPoint(position, R3zero_vector, buf[3]);
  }

  // Close file
  fclose(fp);

  // Return success
  return 1;
}



int R3PointSet::
WriteSDFFile(const char *filename) const
{
 // Open file
  FILE *fp = fopen(filename, "wb");
  if (!fp) {
    RNFail("Unable to open %s\n", filename);
    return 0;
  }

  // Write points
  float buf[4];
  for (int i = 0; i < NPoints(); i++) {
    const R3Point& p = PointPosition(i);
    buf[0] = p[0]; buf[1] = p[1]; buf[2] = p[2];
    buf[3] = PointValue(i);
    if (fwrite(buf, sizeof(float), 4, fp) != (unsigned int) 4) {
      RNFail("Unable to write point %d to %s\n", i, filename);
      return 0;
    }
  }

  // Close file
  fclose(fp);

  // Return success
  return 1;
}



////////////////////////////////////////////////////////////////////////

int R3PointSet::
ReadMeshFile(const char *filename)
{
  // Read mesh file
  R3Mesh mesh;
  if (!mesh.ReadFile(filename)) {
    RNFail("Unable to read %s\n", filename);
    return 0;
  }

  // Create point for each vertex
  for (int i = 0; i < mesh.NVertices(); i++) {
    R3MeshVertex *vertex = mesh.Vertex(i);
    R3Point position = mesh.VertexPosition(vertex);
    R3Vector normal = mesh.VertexNormal(vertex);
    InsertPoint(position, normal);
  }

  // Return success
  return 1;
}



int R3PointSet::
WriteMeshFile(const char *filename) const
{
  // Create mesh
  R3Mesh mesh;
  for (int i = 0; i < NPoints(); i++) {
    mesh.CreateVertex(PointPosition(i), PointNormal(i));
  }

  // Write mesh
  return mesh.WriteFile(filename);
}
