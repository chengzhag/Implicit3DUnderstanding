/* Source file for the GAPS flags class */



/* Include files */

#include "RNBasics.h"



// Namespace

namespace gaps {



/* Public functions */

int 
RNInitFlags()
{
    /* Return success */
    return TRUE;
}



void 
RNStopFlags()
{
}



RNFlags::
RNFlags(void)
    : flags(0)
{
}



} // namespace gaps
