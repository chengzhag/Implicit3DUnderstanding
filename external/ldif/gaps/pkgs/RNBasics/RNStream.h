/* Include file for the GAPS stream class */
#ifndef __RN__STREAM__H__
#define __RN__STREAM__H__



/* Begin namespace */
namespace gaps {



/* Initialization functions */

int RNInitStream();
void RNStopStream();



/* Class definition */

class RNStream {
    public:
        // Constructor functions
        RNStream(void);
        virtual ~RNStream(void);

	// File open/close
        int Open(const char *filename, const char *mode);
        int Close(void);
	
        // Read and write
        size_t Read(void * ptr, size_t size, size_t count);
        size_t Write(void * ptr, size_t size, size_t count);

        // Buffered parsing
        int Scanf(const char *fmt, ...);
        int Printf(const char *fmt, ...);
         
        // Seek and tell
        int Seek(unsigned long long offset, int whence);
        unsigned long long Tell(void);

    private:
        FILE *fp;
        char *buffer_data;
        size_t buffer_size;
        char *filename;
        char *mode;
};



////////////////////////////////////////////////////////////////////////
// Inline utility functions
////////////////////////////////////////////////////////////////////////

inline RNStream *
RNStreamOpen(const char *filename, const char *mode)
{
  // Allocate stream
  RNStream *stream = new RNStream();
  if (!stream) return NULL;

  // Open stream (just like fopen)
  if (!stream->Open(filename, mode)) {
    delete stream;
    return NULL;
  }

  // Return stream
  return stream;
}

  
  
inline int
RNStreamClose(RNStream *stream)
{
  // Close stream (just like fclose)
  return stream->Close();
}

  

inline size_t
RNStreamRead(void * ptr, size_t size, size_t count, RNStream *stream)
{
  // Write to stream (just like fread)
  return stream->Read(ptr, size, count);
}

  

inline size_t
RNStreamWrite(void * ptr, size_t size, size_t count, RNStream *stream)
{
  // Write to stream (just like fwrite)
  return stream->Write(ptr, size, count);
}


  
inline int
RNStreamScanf(RNStream *stream, const char *fmt, ...)
{
  // Read from buffered stream (just like fscanf)
  va_list args;
  va_start(args, fmt);
  int status = stream->Scanf(fmt, args);
  va_end(args);

  // Return status
  return status;
}


  
inline int
RNStreamPrintf(RNStream *stream, const char *fmt, ...)
{
  // Write to buffered stream (just like fprintf)
  va_list args;
  va_start(args, fmt);
  int status = stream->Printf(fmt, args);
  va_end(args);

  // Return status
  return status;
}


  
inline int
RNStreamSeek(RNStream *stream, unsigned long long offset, int whence)
{
  // Seek to position in stream (just like fseek)
  return stream->Seek(offset, whence);
}


  
inline unsigned long long
RNStreamTell(RNStream *stream)
{
  // Tell where current position is in stream (just like ftell)
  return stream->Tell();
}



////////////////////////////////////////////////////////////////////////
// Stream seek constants
////////////////////////////////////////////////////////////////////////

#define RN_STREAM_SEEK_SET SEEK_SET
#define RN_STREAM_SEEK_CUR SEEK_CUR
#define RN_STREAM_SEEK_END SEEK_END



// End namespace
}


// End include guard
#endif
