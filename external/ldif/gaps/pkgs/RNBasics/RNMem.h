/* Include file for GAPS mem utility */
#ifndef __RN__MEMORY__H__
#define __RN__MEMORY__H__



/* Begin namespace */
namespace gaps {



/* Initialization functions */

int RNInitMem();
void RNStopMem();



/* Memory allocation/deallocation functions */

#if FALSE

void *operator new(size_t size);
void *operator new(size_t size, size_t extra);
void operator delete(void *data);

#endif



/* Standard memory manipulation functions */

#define RN_SWAP_BUFFER_SIZE 1024
void RNSwap(void *node1, void *node2, void *buffer, int size);
void RNCopy(const void *src, void *dst, int size);
void RNZero(void *data, int size);
int RNCompare(const void *src1, const void *src2, int size);



/* String copy function */

char *RNStrdup(const char *str);


  
/* Memory usage statistics */

long RNMaxMemoryUsage(void);



// End namespace
}


// End include guard
#endif
