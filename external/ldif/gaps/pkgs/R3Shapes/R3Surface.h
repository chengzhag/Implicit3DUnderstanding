/* Include file for the R3 surface class */
#ifndef __R3__SURFACE__H__
#define __R3__SURFACE__H__



/* Begin namespace */
namespace gaps {



/* Initialization functions */

int R3InitSurface();
void R3StopSurface();



/* Class definition */

class R3Surface : public R3Shape {
    public:
        // Constructors/destructors ???
	R3Surface(void);
	~R3Surface(void);

        // Shape property functions/operators
	const RNBoolean IsSurface(void) const;
};



// End namespace
}


// End include guard
#endif
