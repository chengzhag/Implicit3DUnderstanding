/* Source file for the GAPS stream class */



/* Include files */

#include "RNBasics.h"



// Namespace

namespace gaps {



/* Public functions */

int 
RNInitStream()
{
    /* Return success */
    return TRUE;
}



void 
RNStopStream()
{
}



RNStream::
RNStream(void)
    : fp(NULL),
      buffer_data(NULL),
      buffer_size(0),
      filename(NULL),
      mode(NULL)
{
}



RNStream::
~RNStream(void)
{
    // Close stream
    if (fp) fclose(fp);
  
    // Delete buffer data
    if (buffer_data) free(buffer_data);

    // Delete filename
    if (filename) free(filename);

    // Delete mode
    if (mode) free(mode);

    // Reset everything
    fp = NULL;
    buffer_data = NULL;
    buffer_size = 0;
    filename = NULL;
    mode = NULL;
}




int RNStream::
Open(const char *filename, const char *mode)
{
    // Close previous stream
    if (fp) Close();

    // Copy filename and mode
    if (!filename || !mode) return 0;
    this->filename = strdup(filename);
    if (!this->filename) return 0;
    this->mode = strdup(mode);
    if (!this->mode) { free(this->filename); return 0; }

    // Check if not CNS file
    if (strncmp(filename, "/cns/", 5)) {
        // Open stream for file
        fp = fopen(filename, mode);
        if (!fp) return 0;
    }
    else {
        // Check mode
        if (strchr(mode, '+')) {
            RNAbort("Read/write not supported on /cns files");
        }
        else if (strchr(mode, 'r')) {
            // Get buffer size
            buffer_size = RNFileSize(filename);
            if (buffer_size <= 0) return 0;
          
            // Allocate buffer data
            buffer_data = (char *) malloc(buffer_size);
            if (!buffer_data) {
                buffer_size = 0;
                return 0;
            }
            
            // Read data into buffer

            // Create stream
            fp = fmemopen(buffer_data, buffer_size, mode);
            if (!fp) {
                delete buffer_data;
                buffer_data = NULL;
                buffer_size = 0;
                return 0;
            }
        }
        else if (strchr(mode, 'w')) {
            // Create stream (buffer grows dynamically as write)
            fp = open_memstream(&buffer_data, &buffer_size);
            if (!fp) {
                free(this->filename);
                this->filename = NULL;
                return 0;
            }
        }    
    }

    // Return success
    return 1;
}

  

int RNStream::
Close(void)
{
    // Check stream
    if (!fp) return 0;
  
    // Flush stream
    fflush(fp);

    // Write buffer to file
    if (buffer_data && (buffer_size > 0) &&
        filename && mode && strchr(mode, 'w')) {
    }

    // Close stream
    fclose(fp);

    // Delete buffer data
    if (buffer_data) free(buffer_data);

    // Delete filename
    if (filename) free(filename);

    // Delete mode
    if (mode) free(mode);

    // Reset everything
    fp = NULL;
    buffer_data = NULL;
    buffer_size = 0;
    filename = NULL;
    mode = NULL;

    // Return success
    return 1;
}

  

size_t RNStream::
Read(void *ptr, size_t size, size_t count)
{
    // Check stream
    if (!fp) return 0;

    // Read from stream
    return fread(ptr, size, count, fp);
}


  
size_t RNStream::
Write(void *ptr, size_t size, size_t count)
{
    // Check stream
    if (!fp) return 0;
  
    // Write to stream
    return fwrite(ptr, size, count, fp);
}




int RNStream::
Scanf(const char *fmt, ...)
{
    // Check stream
    if (!fp) return 0;
  
    // Scan from stream
    va_list args;
    va_start(args, fmt);
    int status = vfscanf(fp, fmt, args);
    va_end(args);

    // Return status
    return status;
}



int RNStream::
Printf(const char *fmt, ...)
{
    // Check stream
    if (!fp) return 0;
  
    // Print to stream
    va_list args;
    va_start(args, fmt);
    int status = vfprintf(fp, fmt, args);
    va_end(args);

    // Return status
    return status;
}



int RNStream::
Seek(unsigned long long offset, int whence)
{
    // Check stream
    if (!fp) return 0;

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



unsigned long long RNStream::
Tell(void)
{
    // Check stream
    if (!fp) return 0;

#if (RN_OS == RN_WINDOWS)
    // Windows
    return _ftell64(fp);
#elif (RN_CC_VER == RN_C11)
    // Linux/unix/cygwin etc.
    return ftell(fp);
#else
    // Linux/unix/cygwin etc.
    return ftello(fp);
#endif
}



} // namespace gaps
