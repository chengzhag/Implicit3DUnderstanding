// Include file for file management utilities
#ifndef __RN__FILE__H__
#define __RN__FILE__H__



/* Begin namespace */
namespace gaps {



////////////////////////////////////////////////////////////////////////
// File existence functions
////////////////////////////////////////////////////////////////////////

RNBoolean RNFileExists(const char *filename);



////////////////////////////////////////////////////////////////////////
// File size functions
////////////////////////////////////////////////////////////////////////

unsigned long long RNFileSize(const char *filename);



////////////////////////////////////////////////////////////////////////
// File seek functions
////////////////////////////////////////////////////////////////////////

int RNFileSeek(FILE *fp, unsigned long long offset, int whence);
unsigned long long RNFileTell(FILE *fp);



////////////////////////////////////////////////////////////////////////
// File I/O functions
////////////////////////////////////////////////////////////////////////

int RNReadChar(FILE *fp, char *ptr, int count, int /* swap_endian */);
int RNReadUnsignedChar(FILE *fp, unsigned char *ptr, int count, int /* swap_endian */);
int RNReadShort(FILE *fp, short *ptr, int count, int swap_endian);
int RNReadUnsignedShort(FILE *fp, unsigned short *ptr, int count, int swap_endian);
int RNReadInt(FILE *fp, int *ptr, int count, int swap_endian);
int RNReadUnsignedInt(FILE *fp, unsigned int *ptr, int count, int swap_endian);
int RNReadFloat(FILE *fp, float *ptr, int count, int swap_endian);
int RNReadDouble(FILE *fp, double *ptr, int count, int swap_endian);
int RNReadLongLong(FILE *fp, long long *ptr, int count, int swap_endian);
int RNReadUnsignedLongLong(FILE *fp, unsigned long long *ptr, int count, int swap_endian);
int RNWriteChar(FILE *fp, const char *ptr, int count, int /* swap_endian */);
int RNWriteUnsignedChar(FILE *fp, const unsigned char *ptr, int count, int /* swap_endian */);
int RNWriteShort(FILE *fp, const short *ptr, int count, int swap_endian);
int RNWriteUnsignedShort(FILE *fp, const unsigned short *ptr, int count, int swap_endian);
int RNWriteInt(FILE *fp, const int *ptr, int count, int swap_endian);
int RNWriteUnsignedInt(FILE *fp, const unsigned int *ptr, int count, int swap_endian);
int RNWriteFloat(FILE *fp, const float *ptr, int count, int swap_endian);
int RNWriteDouble(FILE *fp, const double *ptr, int count, int swap_endian);
int RNWriteLongLong(FILE *fp, const long long *ptr, int count, int swap_endian);
int RNWriteUnsignedLongLong(FILE *fp, const unsigned long long *ptr, int count, int swap_endian);



////////////////////////////////////////////////////////////////////////
// Endian byteswapping functions
////////////////////////////////////////////////////////////////////////

void RNSwap2(void *values, int count);
void RNSwap4(void *values, int count);
void RNSwap8(void *values, int count);

  

////////////////////////////////////////////////////////////////////////
// File seek constants
////////////////////////////////////////////////////////////////////////

#define RN_FILE_SEEK_SET SEEK_SET
#define RN_FILE_SEEK_CUR SEEK_CUR
#define RN_FILE_SEEK_END SEEK_END



// End namespace
}


// End include guard
#endif
