// Include file for image class
#ifndef __R2__IMAGE__H__
#define __R2__IMAGE__H__



/* Begin namespace */
namespace gaps {


  
// Class definition

class R2Image {
 public:
  // Constructors
  R2Image(void);
  R2Image(const char *filename);
  R2Image(int width, int height, int ncomponents = 3);
  R2Image(int width, int height, int ncomponents, unsigned char *data);
  R2Image(const R2Grid& red, const R2Grid& green, const R2Grid& blue);
  R2Image(const R2Image& image);
  ~R2Image(void);

  // Accessors
  const unsigned char *Pixels(void) const;
  const unsigned char *Pixels(int row) const;
  const unsigned char *Pixel(int row, int column) const;
  const RNRgb PixelRGB(int row, int column) const;
  int Width(void) const;
  int Height(void) const;
  int Depth(void) const;
  int NComponents(void) const;
  int RowSize(void) const;
  int Size(void) const;

  // Manipulation
  R2Image& operator=(const R2Image& image);
  void Clear(const RNRgb& rgb = RNblack_rgb);
  void Add(const R2Image& image);
  void Subtract(const R2Image& image);
  void Resize(int width, int height, int ncomponents);
  void SetPixel(int row, int column, const unsigned char *pixel);
  void SetPixelRGB(int row, int column, const RNRgb& rgb);

  // File reading/writing
  int ReadFile(const char *filename);
  int ReadBMPFile(const char *filename);
  int ReadPPMFile(const char *filename);
  int ReadPFMFile(const char *filename);
  int ReadJPEGFile(const char *filename);
  int ReadTIFFFile(const char *filename);
  int ReadPNGFile(const char *filename);
  int ReadRAWFile(const char *filename);
  int ReadGRDFile(const char *filename);
  int WriteFile(const char *filename) const;
  int WriteBMPFile(const char *filename) const;
  int WritePPMFile(const char *filename, int ascii = 0) const;
  int WriteRAWFile(const char *filename) const;
  int WriteJPEGFile(const char *filename) const;
  int WriteTIFFFile(const char *filename) const;
  int WritePNGFile(const char *filename) const;

  // Buffer reading/writing
  int ReadBuffer(char *buffer, size_t buffer_length);
  int ReadPNGBuffer(char *buffer, size_t buffer_length);
  int WriteBuffer(char **buffer, size_t *buffer_length) const;
  int WritePNGBuffer(char **buffer, size_t *buffer_length) const;
  
  // Stream reading/writing
  int ReadStream(FILE *fp);
  int ReadPNGStream(FILE *fp);
  int WriteStream(FILE *fp) const;
  int WritePNGStream(FILE *fp) const;
  
  // Capture functions
  void Capture(void);

  // Draw functions
  void Draw(int x = 0, int y = 0) const;

public:
  // For backward compatibility
  int Read(const char *filename) { return ReadFile(filename); };
  int ReadBMP(const char *filename) { return ReadBMPFile(filename); };
  int ReadPPM(const char *filename) { return ReadPPMFile(filename); };
  int ReadPFM(const char *filename) { return ReadPFMFile(filename); };
  int ReadJPEG(const char *filename) { return ReadJPEGFile(filename); };
  int ReadTIFF(const char *filename) { return ReadTIFFFile(filename); };
  int ReadPNG(const char *filename) { return ReadPNGFile(filename); };
  int ReadRAW(const char *filename) { return ReadRAWFile(filename); };
  int ReadGRD(const char *filename) { return ReadGRDFile(filename); };
  int Write(const char *filename) const { return WriteFile(filename); };
  int WriteBMP(const char *filename) const { return WriteBMPFile(filename); };
  int WritePPM(const char *filename, int ascii = 0) const { return WritePPMFile(filename, ascii); };
  int WriteRAW(const char *filename) const { return WriteRAWFile(filename); };
  int WriteJPEG(const char *filename) const { return WriteJPEGFile(filename); };
  int WriteTIFF(const char *filename) const { return WriteTIFFFile(filename); };
  int WritePNG(const char *filename) const { return WritePNGFile(filename); };

private:
  int width;
  int height;
  int ncomponents;
  int rowsize;
  unsigned char *pixels;
};



// Inline functions

inline int R2Image::
Width(void) const
{
  // Return width
  return width;
}



inline int R2Image::
Height(void) const
{
  // Return height
  return height;
}



inline int R2Image::
Depth(void) const
{
  // Return number of bytes per pixel
  return ncomponents;
}



inline int R2Image::
NComponents(void) const
{
  // Return number of bytes per pixel
  return Depth();
}



inline int R2Image::
RowSize(void) const
{
  // Return size of row in bytes
  return rowsize;
}



inline int R2Image::
Size(void) const
{
  // Return size of image in bytes
  return rowsize * height;
}



inline const unsigned char *R2Image::
Pixels(void) const
{
  // Return pixels pointer (pixels start at lower-left)
  return pixels;
}



inline const unsigned char *R2Image::
Pixels(int y) const
{
  // Return pixels pointer for row y
  return &pixels[y*rowsize];
}



inline const unsigned char *R2Image::
Pixel(int x, int y) const
{
  // Return pixel value at (x,y)
  return &pixels[y*rowsize + x*ncomponents];
}


  
inline int R2Image::
ReadBuffer(char *buffer, size_t buffer_length)
{
  // Read from PNG buffer
  return ReadPNGBuffer(buffer, buffer_length);
}


  
inline int R2Image::
WriteBuffer(char **buffer, size_t *buffer_length) const
{
  // Write to PNG buffer
  return WritePNGBuffer(buffer, buffer_length);
}



inline int R2Image::
ReadStream(FILE *fp)
{
  // Read from PNG stream
  return ReadPNGStream(fp);
}

  
inline int R2Image::
WriteStream(FILE *fp) const
{
  // Write to PNG stream
  return WritePNGStream(fp);
}

  

  
// End namespace
}


// End include guard
#endif
