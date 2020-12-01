/* Include file for R2 shapes module */
#ifndef __R2__SHAPES__H__
#define __R2__SHAPES__H__



/* Dependency include files */

#include "RNBasics/RNBasics.h"



/* Class declarations */

namespace gaps {
class R2Vector;
class R2Point;
class R2Line;
class R2Ray;
class R2Span;
class R2Halfspace;
class R3Matrix;
class R2Diad;
class R2CoordSystem;
class R2Transformation;
class R2Affine;
class R2Arc;
class R2Polyline;
class R2Box;
class R2Circle;
class R2Polygon;
class R2Grid;
class R2PixelDatabase;
}


/* Image include files */

#include "R2Image.h"



/* Primitive shape include files */

#include "R2Vector.h"
#include "R2Point.h"
#include "R2Line.h"
#include "R2Ray.h"
#include "R2Span.h"
#include "R2Halfspace.h"



/* Transformation include files */

#include "R3Matrix.h"
#include "R2Diad.h"
#include "R2Crdsys.h"
#include "R2Xform.h"
#include "R2Affine.h"



/* Abstract shape include file */

#include "R2Shape.h"



/* Solid shapes include files */

#include "R2Solid.h"
#include "R2Box.h"
#include "R2Circle.h"
#include "R2Polygon.h"



/* Curve shapes include files */

#include "R2Curve.h"
#include "R2Arc.h"
#include "R2Polyline.h"



/* Image/grid include files */

#include "R2Grid.h"
#include "R2PixelDatabase.h"


/* Shape relationship include files */

#include "R2Perp.h"
#include "R2Parall.h"
#include "R2Dist.h"
#include "R2Cont.h"
#include "R2Isect.h"
#include "R2Relate.h"
#include "R2Align.h"



/* Closest point search include files */

#include "R2Kdtree.h"



/* Shape utility include files */

#include "R2Draw.h"
#include "R2Io.h"



/* Initialization functions */

namespace gaps {
int R2InitShapes(void);
void R2StopShapes(void);
}



#endif







