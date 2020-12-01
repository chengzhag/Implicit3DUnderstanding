/* Include file for input/output utility */
#ifndef __R2__IO__H__
#define __R2__IO__H__



/* Begin namespace */
namespace gaps {


  
/* Function declarations */

RNArray<R2Span *> *R2ReadSpans(RNArray<R2Span *>& spans, const char *filename);
RNArray<R2Span *> *R2ReadWiseFile(RNArray<R2Span *>& spans, const char *filename);
RNArray<R2Span *> *R2ReadXFigFile(RNArray<R2Span *>& spans, const char *filename);



// End namespace
}


// End include guard
#endif
