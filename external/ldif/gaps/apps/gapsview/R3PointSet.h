// Include file for 3D point set


// Include files

#include <vector>



// Class definition

class R3PointSet {
public:
  // Constructor/deconstructor
  R3PointSet(void);
  ~R3PointSet(void);

  // Point access
  int NPoints(void) const;
  const R3Point& PointPosition(int index) const;
  const R3Vector& PointNormal(int index) const;
  RNScalar PointValue(int index) const;

  // Point set properties
  R3Point Centroid(void) const;
  const R3Box& BBox(void) const;

  // Point queries
  int FindPointIndex(const R3Point& position, double tolerance = 0) const;

  // Point manipulation
  void RemovePoint(int index);
  void InsertPoint(const R3Point& position, const R3Vector& normal, RNScalar value = 0);
  void SetPointPosition(int index, const R3Point& position);
  void SetPointNormal(int index, const R3Vector& normal);
  void SetPointValue(int index, RNScalar value);

  // Input/output
  int ReadFile(const char *filename);
  int ReadASCIIFile(const char *filename);
  int ReadBinaryFile(const char *filename);
  int ReadMeshFile(const char *filename);
  int ReadSDFFile(const char *filename);
  int WriteFile(const char *filename) const;
  int WriteASCIIFile(const char *filename) const;
  int WriteBinaryFile(const char *filename) const;
  int WriteSDFFile(const char *filename) const;
  int WriteMeshFile(const char *filename) const;

private:
  std::vector<R3Point> points;
  std::vector<R3Vector> normals;
  std::vector<RNScalar> values;
  R3Box bbox;
};



// Inline functions

inline int R3PointSet::
NPoints(void) const
{
  // Return number of points
  return (int) points.size();
}



inline const R3Point& R3PointSet::
PointPosition(int index) const
{
  // Return position of point with given index
  return points[index];
}



inline const R3Vector& R3PointSet::
PointNormal(int index) const
{
  // Return normal of point with given index
  return normals[index];
}



inline RNScalar R3PointSet::
PointValue(int index) const
{
  // Return value of point with given index
  return values[index];
}



inline const R3Box& R3PointSet::
BBox(void) const
{
  // Return bounding box of point set
  return bbox;
}


