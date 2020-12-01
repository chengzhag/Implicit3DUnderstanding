// Source file for file management utilities


// Include files
#include "RNBasics.h"



// Namespace

namespace gaps {



////////////////////////////////////////////////////////////////////////
// FILE EXISTENCE FUNCTIONS
////////////////////////////////////////////////////////////////////////

RNBoolean
RNFileExists(const char *filename)
{
    // Return whether or not file exists (and is readable)
    FILE *fp = fopen(filename, "rb");
    if (!fp) return FALSE;
    fclose(fp); 
    return TRUE;
}



////////////////////////////////////////////////////////////////////////
// FILE SIZE FUNCTIONS
////////////////////////////////////////////////////////////////////////

unsigned long long
RNFileSize(const char *filename)
{
#if 1
    // Return size of file in bytes
    FILE *fp = fopen(filename, "r");
    if (!fp) return 0;
    RNFileSeek(fp, RN_FILE_SEEK_END, 0);
    unsigned long long file_size = RNFileTell(fp);
    fclose(fp);
    return file_size;
#else
    // Requires "unistd.h"
    struct stat stat_buf;
    int rc = fstat(filename, &stat_buf);
    if (rc == 0) return 0;
    else return stat_buf.st_size;
#endif
}



////////////////////////////////////////////////////////////////////////
// FILE SEEK FUNCTIONS
////////////////////////////////////////////////////////////////////////

int 
RNFileSeek(FILE *fp, unsigned long long offset, int whence)
{
#if (RN_OS == RN_WINDOWS)
    // Windows
    if (_fseek64(fp, offset, whence) == 0) return 1;
    else return 0;
#elif (RN_CC_VER == RN_C11)
    // Linux/unix/cygwin etc.
    if (fseek(fp, offset, whence) == 0) return 1;
    else return 0;
#else
    // Linux/unix/cygwin etc.
    if (fseeko(fp, offset, whence) == 0) return 1;
    else return 0;
#endif
}



unsigned long long 
RNFileTell(FILE *fp)
{
#if (RN_OS == RN_WINDOWS)
    return _ftell64(fp);
#elif (RN_CC_VER == RN_C11)
    return ftell(fp);
#else
    // Linux/unix/cygwin etc.
    return ftello(fp);
#endif
}
  


////////////////////////////////////////////////////////////////////////
// FILE I/O FUNCTIONS
////////////////////////////////////////////////////////////////////////

int 
RNReadChar(FILE *fp, char *ptr, int count, int /* swap_endian */)
{
  // Read the values
  if (fread(ptr, sizeof(char), count, fp) != (size_t) count) {
    RNFail("Unable to read char from file\n");
    return 0;
  }

  // Return success
  return 1;
}



int 
RNReadUnsignedChar(FILE *fp, unsigned char *ptr, int count, int /* swap_endian */)
{
  // Read the values
  if (fread(ptr, sizeof(unsigned char), count, fp) != (size_t) count) {
    RNFail("Unable to read unsigned char from file\n");
    return 0;
  }

  // Return success
  return 1;
}



int 
RNReadShort(FILE *fp, short *ptr, int count, int swap_endian)
{
  // Read the values
  if (fread(ptr, sizeof(short), count, fp) != (size_t) count) {
    RNFail("Unable to read short from file\n");
    return 0;
  }

  // Swap endian
  if (swap_endian) RNSwap2(ptr, count);

  // Return success
  return 1;
}



int 
RNReadUnsignedShort(FILE *fp, unsigned short *ptr, int count, int swap_endian)
{
  // Read the values
  if (fread(ptr, sizeof(unsigned short), count, fp) != (size_t) count) {
    RNFail("Unable to read unsigned short from file\n");
    return 0;
  }

  // Swap endian
  if (swap_endian) RNSwap2(ptr, count);

  // Return success
  return 1;
}



int 
RNReadInt(FILE *fp, int *ptr, int count, int swap_endian)
{
  // Read the values
  if (fread(ptr, sizeof(int), count, fp) != (size_t) count) {
    RNFail("Unable to read integer from file\n");
    return 0;
  }

  // Swap endian
  if (swap_endian) RNSwap4(ptr, count);

  // Return success
  return 1;
}



int 
RNReadUnsignedInt(FILE *fp, unsigned int *ptr, int count, int swap_endian)
{
  // Read the values
  if (fread(ptr, sizeof(unsigned int), count, fp) != (size_t) count) {
    RNFail("Unable to read unsigned integer from file\n");
    return 0;
  }

  // Swap endian
  if (swap_endian) RNSwap4(ptr, count);

  // Return success
  return 1;
}



int
RNReadFloat(FILE *fp, float *ptr, int count, int swap_endian)
{
  // Read the values
  if (fread(ptr, sizeof(float), count, fp) != (size_t) count) {
    RNFail("Unable to read float from file\n");
    return 0;
  }

  // Swap endian
  if (swap_endian) RNSwap4(ptr, count);

  // Return success
  return 1;
}



int
RNReadDouble(FILE *fp, double *ptr, int count, int swap_endian)
{
  // Read the values
  if (fread(ptr, sizeof(double), count, fp) != (size_t) count) {
    RNFail("Unable to read double from file\n");
    return 0;
  }

  // Swap endian
  if (swap_endian) RNSwap8(ptr, count);

  // Return success
  return 1;
}



int
RNReadLongLong(FILE *fp, long long *ptr, int count, int swap_endian)
{
  // Read the values
  if (fread(ptr, sizeof(long long), count, fp) != (size_t) count) {
    RNFail("Unable to read long long from file\n");
    return 0;
  }

  // Swap endian
  if (swap_endian) RNSwap8(ptr, count);

  // Return success
  return 1;
}



int
RNReadUnsignedLongLong(FILE *fp, unsigned long long *ptr, int count, int swap_endian)
{
  // Read the values
  if (fread(ptr, sizeof(unsigned long long), count, fp) != (size_t) count) {
    RNFail("Unable to read unsigned long long from file\n");
    return 0;
  }

  // Swap endian
  if (swap_endian) RNSwap8(ptr, count);

  // Return success
  return 1;
}



int
RNWriteChar(FILE *fp, const char *ptr, int count, int /* swap_endian */)
{
  // Write the values
  if (fwrite(ptr, sizeof(char), count, fp) != (size_t) count) {
    RNFail("Unable to write integer to file\n");
    return 0;
  }

  // Return success
  return 1;
}



int
RNWriteUnsignedChar(FILE *fp, const unsigned char *ptr, int count, int /* swap_endian */)
{
  // Write the values
  if (fwrite(ptr, sizeof(unsigned char), count, fp) != (size_t) count) {
    RNFail("Unable to write unsigned char to file\n");
    return 0;
  }

  // Return success
  return 1;
}



int
RNWriteShort(FILE *fp, const short *ptr, int count, int swap_endian)
{
  // Swap endian
  if (swap_endian) RNSwap2((void *) ptr, count);

  // Write the values
  int status = 1;
  if (fwrite(ptr, sizeof(short), count, fp) != (size_t) count) {
    RNFail("Unable to write short to file\n");
    status = 0;
  }

  // Swap endian back
  if (swap_endian) RNSwap2((void *) ptr, count);

  // Return status
  return status;
}



int
RNWriteUnsignedShort(FILE *fp, const unsigned short *ptr, int count, int swap_endian)
{
  // Swap endian
  if (swap_endian) RNSwap2((void *) ptr, count);

  // Write the values
  int status = 1;
  if (fwrite(ptr, sizeof(unsigned short), count, fp) != (size_t) count) {
    RNFail("Unable to write unsigned short to file\n");
    status = 0;
  }

  // Swap endian back
  if (swap_endian) RNSwap2((void *) ptr, count);

  // Return status
  return status;
}



int
RNWriteInt(FILE *fp, const int *ptr, int count, int swap_endian)
{
  // Swap endian
  if (swap_endian) RNSwap4((void *) ptr, count);

  // Write the values
  int status = 1;
  if (fwrite(ptr, sizeof(int), count, fp) != (size_t) count) {
    RNFail("Unable to write integer to file\n");
    status = 0;
  }

  // Swap endian back
  if (swap_endian) RNSwap4((void *) ptr, count);

  // Return status
  return status;
}



int
RNWriteUnsignedInt(FILE *fp, const unsigned int *ptr, int count, int swap_endian)
{
  // Swap endian
  if (swap_endian) RNSwap4((void *) ptr, count);

  // Write the values
  int status = 1;
  if (fwrite(ptr, sizeof(unsigned int), count, fp) != (size_t) count) {
    RNFail("Unable to write unsigned integer to file\n");
    status = 0;
  }

  // Swap endian back
  if (swap_endian) RNSwap4((void *) ptr, count);

  // Return status
  return status;
}



int
RNWriteFloat(FILE *fp, const float *ptr, int count, int swap_endian)
{
  // Swap endian
  if (swap_endian) RNSwap4((void *) ptr, count);

  // Write the values
  int status = 1;
  if (fwrite(ptr, sizeof(float), count, fp) != (size_t) count) {
    RNFail("Unable to write float to file\n");
    status = 0;
  }

  // Swap endian back
  if (swap_endian) RNSwap4((void *) ptr, count);

  // Return status
  return status;
}



int
RNWriteDouble(FILE *fp, const double *ptr, int count, int swap_endian)
{
  // Swap endian
  if (swap_endian) RNSwap8((void *) ptr, count);

  // Write the values
  int status = 1;
  if (fwrite(ptr, sizeof(double), count, fp) != (size_t) count) {
    RNFail("Unable to write double to file\n");
    status = 0;
  }

  // Swap endian back
  if (swap_endian) RNSwap8((void *) ptr, count);

  // Return status
  return status;
}



int
RNWriteLongLong(FILE *fp, const long long *ptr, int count, int swap_endian)
{
  // Swap endian
  if (swap_endian) RNSwap8((void *) ptr, count);

  // Write the values
  int status = 1;
  if (fwrite(ptr, sizeof(long long), count, fp) != (size_t) count) {
    RNFail("Unable to write long long to file\n");
    status = 0;
  }

  // Swap endian back
  if (swap_endian) RNSwap8((void *) ptr, count);

  // Return status
  return status;
}



int
RNWriteUnsignedLongLong(FILE *fp, const unsigned long long *ptr, int count, int swap_endian)
{
  // Swap endian
  if (swap_endian) RNSwap8((void *) ptr, count);

  // Write the values
  int status = 1;
  if (fwrite(ptr, sizeof(unsigned long long), count, fp) != (size_t) count) {
    RNFail("Unable to write unsigned long long to file\n");
    status = 0;
  }

  // Swap endian back
  if (swap_endian) RNSwap8((void *) ptr, count);

  // Return status
  return status;
}



////////////////////////////////////////////////////////////////////////
// ENDIAN BYTE-SWAPPING FUNCTIONS
////////////////////////////////////////////////////////////////////////

void
RNSwap2(void *values, int count)
{
  // Swap endian of 2-byte data type
  unsigned short *y = (unsigned short *) values;
  for (int i = 0; i < count; i++) {
    unsigned short x = y[i];
    y[i] = (x<<8) | (x>>8);
  }
}



void
RNSwap4(void *values, int count)
{
  // Swap endian of 4-byte data type
  unsigned int *y = (unsigned int *) values;
  for (int i = 0; i < count; i++) {
    unsigned int x = y[i];
    y[i] = 
       (x<<24) | 
      ((x<<8) & 0x00FF0000) | 
      ((x>>8) & 0x0000FF00) | 
       (x>>24);
  }
}



void
RNSwap8(void *values, int count)
{
  // Swap endian
  unsigned long long *y = (unsigned long long *) values;
  for (int i = 0; i < count; i++) {
    unsigned long long x = y[i];
    y[i] = 
       (x<<56) | 
      ((x<<40) & 0x00FF000000000000ULL) |
      ((x<<24) & 0x0000FF0000000000ULL) |
      ((x<<8)  & 0x000000FF00000000ULL) |
      ((x>>8)  & 0x00000000FF000000ULL) |
      ((x>>24) & 0x0000000000FF0000ULL) |
      ((x>>40) & 0x000000000000FF00ULL) |
       (x>>56);
  }
}



} // namespace gaps
