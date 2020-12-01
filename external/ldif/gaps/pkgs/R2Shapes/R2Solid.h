/* Include file for the R2 solid class */
#ifndef __R2__SOLID__H__
#define __R2__SOLID__H__



/* Begin namespace */
namespace gaps {



/* Initialization functions */

int R2InitSolid();
void R2StopSolid();



/* Class definition */

class R2Solid : public R2Shape {
    public:
        // Constructors/destructors ???
	R2Solid(void);
	~R2Solid(void);

        // Shape property functions/operators
	const RNBoolean IsSolid(void) const;
};



// End namespace
}


// End include guard
#endif
