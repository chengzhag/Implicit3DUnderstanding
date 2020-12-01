/* Include file for the R3 solid class */
#ifndef __R3__SOLID__H__
#define __R3__SOLID__H__



/* Begin namespace */
namespace gaps {



/* Initialization functions */

int R3InitSolid();
void R3StopSolid();



/* Class definition */

class R3Solid : public R3Shape {
    public:
        // Constructors/destructors ???
	R3Solid(void);
	~R3Solid(void);

        // Shape property functions/operators
	const RNBoolean IsSolid(void) const;
};



// End namespace
}


// End include guard
#endif
