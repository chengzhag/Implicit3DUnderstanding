// Header file for GAPS scalar pixel database class
#ifndef __R2__PIXEL__DATABASE__H__
#define __R2__PIXEL__DATABASE__H__



/* Begin namespace */
namespace gaps {


  
// Class definition

class R2PixelDatabase {
public:
  // Constructors
  R2PixelDatabase(void);
  R2PixelDatabase(const R2PixelDatabase& database);
  virtual ~R2PixelDatabase(void);

  // Property access functions
  const char *Name(void) const;
  RNBoolean IsOpen(void) const;
  
  // Entry access functions
  int NKeys(void) const;
  const char *Key(int k) const;
  virtual int FindImage(const char *key, R2Image *image) const;
  virtual int FindGrid(const char *key, R2Grid *grid) const;

  // Property manipulation functions
  virtual void SetName(const char *name);

  // Entry manipulation functions
  virtual int InsertImage(const char *key, const R2Image& image);
  virtual int InsertGrid(const char *key, const R2Grid& grid,
    double offset = 0, double scale = 1, double exponent = 1,
    RNBoolean apply_transformation = TRUE);
  virtual int Remove(const char *key);

  // I/O functions
  virtual int OpenFile(const char *filename, const char *rwaccess = NULL);
  virtual int CloseFile(void);

public:
  // For backward compatibility
  int NEntries() const;

protected:
  // Internal I/O functions
  virtual int ReadHeader(FILE *fp);
  virtual int WriteHeader(FILE *fp, int swap_endian);
  virtual int ReadEntries(FILE *fp, int swap_endian);
  virtual int WriteEntries(FILE *fp, int swap_endian);
  
private:
  FILE *fp;
  char *name;
  char *filename;
  char *rwaccess;
  unsigned int major_version;
  unsigned int minor_version;
  unsigned int swap_endian;
  unsigned int entries_count;
  unsigned long long entries_seek;
  RNSymbolTable<struct R2PixelDatabaseEntry> map;
  std::vector<std::string> keys;
};



////////////////////////////////////////////////////////////////////////
// INLINE FUNCTION DEFINITIONS
////////////////////////////////////////////////////////////////////////

inline const char *R2PixelDatabase::
Name(void) const
{
  // Return name
  return name;
}


  
inline int R2PixelDatabase::
NKeys(void) const
{
  // Return number of entries
  assert(entries_count == keys.size());
  assert((int) entries_count == map.NEntries());
  return (int) entries_count;
}


  
inline const char *R2PixelDatabase::
Key(int k) const
{
  // Return kth key
  return keys[k].c_str();
}
  
  
  
inline RNBoolean R2PixelDatabase::
IsOpen(void) const
{
  // Return whether database is open
  return (fp) ? TRUE : FALSE;
}



inline int R2PixelDatabase::
NEntries(void) const
{
  // Return number of entries
  // Same as NKeys
  // Remains for backward compatibility
  return NKeys();
}


  
// End namespace
}


// End include guard
#endif
