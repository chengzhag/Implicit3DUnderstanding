// Include file for mesh search tree class
#ifndef __R3__MESH__SEARCH__TREE__H__
#define __R3__MESH__SEARCH__TREE__H__



/* Begin namespace */
namespace gaps {



// Node declaration

class R3MeshSearchTreeFace;
class R3MeshSearchTreeNode;



// Class declaration

class R3MeshSearchTree {
public:
  // Constructor/destructors
  R3MeshSearchTree(R3Mesh *mesh);
  ~R3MeshSearchTree(void);

  // Property functions
  R3Mesh *Mesh(void) const;
  const R3Box& BBox(void) const;

  // Insert/delete functions
  void InsertFace(R3MeshFace *face);
  void Empty(void);

  // Find mesh feature closest to a query point
  void FindClosest(const R3Point& query, R3MeshIntersection& closest,
    RNScalar min_distance = 0, RNScalar max_distance = RN_INFINITY,
    int (*IsCompatible)(const R3Point&, const R3Vector&, R3Mesh *, R3MeshFace *, void *) = NULL, 
    void *compatible_data = NULL);

  // Find mesh feature closest to a query point and normal
  void FindClosest(const R3Point& query, const R3Vector& normal, R3MeshIntersection& closest,
    RNScalar min_distance = 0, RNScalar max_distance = RN_INFINITY, 
    int (*IsCompatible)(const R3Point&, const R3Vector&, R3Mesh *, R3MeshFace *, void *) = NULL, 
    void *compatible_data = NULL);

  // Find all mesh features with distance from a query point
  void FindAll(const R3Point& query, RNArray<R3MeshIntersection *>& hits,
    RNScalar min_distance = 0, RNScalar max_distance = RN_INFINITY,
    int (*IsCompatible)(const R3Point&, const R3Vector&, R3Mesh *, R3MeshFace *, void *) = NULL, 
    void *compatible_data = NULL);

  // Find all mesh features with distance from a query point (and with a compatible normal)
  void FindAll(const R3Point& query, const R3Vector& normal, RNArray<R3MeshIntersection *>& hits,
    RNScalar min_distance = 0, RNScalar max_distance = RN_INFINITY,
    int (*IsCompatible)(const R3Point&, const R3Vector&, R3Mesh *, R3MeshFace *, void *) = NULL, 
    void *compatible_data = NULL);

  // Find first ray intersection
  void FindIntersection(const R3Ray& ray, R3MeshIntersection& closest,
    RNScalar min_t = 0, RNScalar max_t = RN_INFINITY,
    int (*IsCompatible)(const R3Point&, const R3Vector&, R3Mesh *, R3MeshFace *, void *) = NULL, 
    void *compatible_data = NULL);

  // Find all mesh faces intersecting shape
  void FindAll(const R3Shape& shape, RNArray<R3MeshIntersection *>& hits);

  // Visualization/debugging functions
  int NNodes(void) const;
  void Outline(void) const;
  void Print(void) const;

public:
  // Internal manipulations functions
  void Empty(R3MeshSearchTreeNode *node);

  // Internal insert functions
  void InsertFace(R3MeshSearchTreeFace *face, R3MeshSearchTreeNode *node, const R3Box& node_box, int depth);

  // Internal closest point search functions
  void FindClosest(const R3Point& query, const R3Vector& normal, R3MeshIntersection& closest, 
    RNScalar min_distance_squared, RNScalar& max_distance_squared, 
    int (*IsCompatible)(const R3Point&, const R3Vector&, R3Mesh *, R3MeshFace *, void *), void *compatible_data,
    R3MeshSearchTreeNode *node, const R3Box& node_box) const;
  void FindClosest(const R3Point& query, const R3Vector& normal, R3MeshIntersection& closest, 
    RNScalar min_distance_squared, RNScalar& max_distance_squared, 
    int (*IsCompatible)(const R3Point&, const R3Vector&, R3Mesh *, R3MeshFace *, void *), void *compatible_data,
    R3MeshFace *face) const;

  // Internal all point search functions
  void FindAll(const R3Point& query, const R3Vector& normal, RNArray<R3MeshIntersection *>& hits,
    RNScalar min_distance_squared, RNScalar max_distance_squared, 
    int (*IsCompatible)(const R3Point&, const R3Vector&, R3Mesh *, R3MeshFace *, void *), void *compatible_data,
    R3MeshSearchTreeNode *node, const R3Box& node_box) const;
  void FindAll(const R3Point& query, const R3Vector& normal, RNArray<R3MeshIntersection *>& hits,
    RNScalar min_distance_squared, RNScalar max_distance_squared, 
    int (*IsCompatible)(const R3Point&, const R3Vector&, R3Mesh *, R3MeshFace *, void *), void *compatible_data,
    R3MeshFace *face) const;

  // Internal ray intersection search functions
  void FindIntersection(const R3Ray& ray, R3MeshIntersection& closest, 
    RNScalar min_t, RNScalar& max_t, 
    int (*IsCompatible)(const R3Point&, const R3Vector&, R3Mesh *, R3MeshFace *, void *), void *compatible_data,
    R3MeshSearchTreeNode *node, const R3Box& node_box) const;
  void FindIntersection(const R3Ray& ray, R3MeshIntersection& closest, 
    RNScalar min_t, RNScalar& max_t, 
    int (*IsCompatible)(const R3Point&, const R3Vector&, R3Mesh *, R3MeshFace *, void *), void *compatible_data,
    R3MeshFace *face) const;

  // Internal shape intersection functions
  void FindAll(const R3Shape& shape, RNArray<R3MeshIntersection *>& hits,
    R3MeshSearchTreeNode *node, const R3Box& node_box) const;

  // Internal visualization and debugging functions
  void Outline(R3MeshSearchTreeNode *node, const R3Box& node_box) const;
  int Print(R3MeshSearchTreeNode *node, int depth) const;

  // Internal utility functions
  RNScalar DistanceSquared(const R3Point& query, const R3Box& box, RNScalar max_distance_squared) const;
  RNScalar DistanceSquared(const R3Point& query, const R3Point& point) const;

  // Not implemented
  R3MeshSearchTree(const R3MeshSearchTree& tree);
  R3MeshSearchTree& operator=(const R3MeshSearchTree& tree);
  
public:
  // Internal data
  R3Mesh *mesh;
  R3MeshSearchTreeNode *root;
  int nnodes;
  RNMark mark;
};



////////////////////////////////////////////////////////////////////////
// Distance functions
////////////////////////////////////////////////////////////////////////

// Distance functions
RNLength R3Distance(const R3MeshSearchTree& tree1, const R3MeshSearchTree& tree2, 
  RNScalar min_distance, RNScalar max_distance,
  int (*IsCompatible)(R3Mesh *, R3MeshFace *, R3Mesh *, R3MeshFace *, void *) = NULL, void *compatible_data = NULL,
  R3MeshIntersection *closest1 = NULL, R3MeshIntersection *closest2 = NULL); 



////////////////////////////////////////////////////////////////////////
// Inline functions
////////////////////////////////////////////////////////////////////////

inline R3Mesh *R3MeshSearchTree::
Mesh(void) const
{
  // Return mesh
  return mesh;
}



inline const R3Box& R3MeshSearchTree::
BBox(void) const
{
  // Return bounding box of the whole KD tree
  return mesh->BBox();
}



// End namespace
}


// End include guard
#endif
