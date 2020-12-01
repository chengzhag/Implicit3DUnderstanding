/* Include file for external stuff */
#ifndef __RN__EXTERN__H__
#define __RN__EXTERN__H__



/* Usual C include files */

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <assert.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <limits.h>



/* Standard library include files */

#include <string>
#include <map>



/* Machine dependent include files */

#if (RN_OS == RN_WINDOWS)
#   include <float.h>
#   include <windows.h>
#else 
#   include <float.h>
#   include <sys/time.h>
#   include <sys/resource.h>
#endif



// End include guard
#endif
