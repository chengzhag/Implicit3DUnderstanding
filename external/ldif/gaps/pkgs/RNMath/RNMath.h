// Include file for RNMath package
#ifndef __RN__MATH__H__
#define __RN__MATH__H__



// Dependency include files

#include "RNBasics/RNBasics.h"



// Declarations

namespace gaps {
class RNVector;
class RNMatrix;
class RNDenseMatrix;
class RNPolynomial;
class RNPolynomialTerm;
class RNAlgebraic;
typedef RNAlgebraic RNExpression;
class RNEquation;
class RNSystemOfEquations;
}



// Matrix classes

#include "RNLapack.h"
#include "RNVector.h"
#include "RNMatrix.h"
#include "RNDenseMatrix.h"
#include "RNDenseLUMatrix.h"


// Expression and equation classes

#include "RNPolynomial.h"
#include "RNAlgebraic.h"
#include "RNEquation.h"
#include "RNSystemOfEquations.h"



#endif
