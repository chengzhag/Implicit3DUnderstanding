/* Include file for R2 draw utility */
#ifndef __R2__DRAW__H__
#define __R2__DRAW__H__



/* Begin namespace */
namespace gaps {


  
inline void 
R2LoadPoint(const R2Point& point)
{
    // Load vertex 
    R2LoadPoint(point.Coords());
}



// End namespace
}


// End include guard
#endif
