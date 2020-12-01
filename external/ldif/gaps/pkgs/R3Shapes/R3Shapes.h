/* Include file for R3 shapes module */
#ifndef __R3__SHAPES__H__
#define __R3__SHAPES__H__



/* Dependency include files */

#include "R2Shapes/R2Shapes.h"



/* Class declarations */

namespace gaps {
class R3Vector;
class R3Point;
class R3Line;
class R3Ray;
class R3Span;
class R3Plane;
class R3Halfspace;
class R4Matrix;
class R3Quaternion;
class R3Triad;
class R3CoordSystem;
class R3Transformation;
class R3Affine;
class R3Shape;
class R3Solid;
class R3Box;
class R3OrientedBox;
class R3Cylinder;
class R3Cone;
class R3Sphere;
class R3Ellipsoid;
class R3Surface;
class R3Triangle;
class R3TriangleArray;
class R3Circle;
class R3Ellipse;
class R3Rectangle;
class R3Mesh;
class R3Curve;
class R3Polyline;
class R3CatmullRomSpline;
class R3PlanarGrid;
class R3Grid;
}



/* Geometry basics include files */

#include "R3Base.h"



/* Primitive include files */

#include "R3Vector.h"
#include "R3Point.h"
#include "R3Line.h"
#include "R3Ray.h"
#include "R3Span.h"
#include "R3Plane.h"
#include "R3Halfspace.h"



/* Transformation include files */

#include "R4Matrix.h"
#include "R3Quaternion.h"
#include "R3Triad.h"
#include "R3Crdsys.h"
#include "R3Xform.h"
#include "R3Affine.h"



/* Abstract shape include files */

#include "R3Shape.h"



/* Some solid shapes include files */

#include "R3Solid.h"
#include "R3Box.h"        



/* Surface shapes include files */

#include "R3Surface.h"
#include "R3Triangle.h"
#include "R3TriangleArray.h"
#include "R3Circle.h"
#include "R3Ellipse.h"
#include "R3Rectangle.h"
#include "R3Mesh.h"
#include "R3PlanarGrid.h"        



/* Surface shapes include files */

#include "R3Curve.h"
#include "R3Polyline.h"
#include "R3CatmullRomSpline.h"



/* More solid shapes include files */

#include "R3OrientedBox.h"        
#include "R3Cylinder.h"
#include "R3Cone.h"
#include "R3Sphere.h"
#include "R3Ellipsoid.h"
#include "R3Grid.h"        



/* Shape relationship include files */

#include "R3Perp.h"
#include "R3Parall.h"
#include "R3Dist.h"
#include "R3Cont.h"
#include "R3Isect.h"
#include "R3Relate.h"
#include "R3Align.h"
#include "R3Kdtree.h"



/* Mesh utility include files */

#include "R3MeshSearchTree.h"
#include "R3MeshProperty.h"
#include "R3MeshPropertySet.h"



/* Shape utility include files */

#include "R3Draw.h"



/* Initialization functions */

namespace gaps {
int R3InitShapes(void);
void R3StopShapes(void);
}



#endif








