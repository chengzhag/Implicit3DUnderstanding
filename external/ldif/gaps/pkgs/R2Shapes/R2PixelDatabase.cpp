/* Source file for the R2 pixel database class */



////////////////////////////////////////////////////////////////////////
// INCLUDE FILES
////////////////////////////////////////////////////////////////////////

#include "R2Shapes.h"
#include <algorithm>



////////////////////////////////////////////////////////////////////////
// Namespace
////////////////////////////////////////////////////////////////////////

namespace gaps {



////////////////////////////////////////////////////////////////////////
// Internal variables
////////////////////////////////////////////////////////////////////////

static const unsigned int current_major_version = 0;
static const unsigned int current_minor_version = 1;

static const unsigned int R2_PIXEL_DATABASE_3_8_PNG_FORMAT = 0;
static const unsigned int R2_PIXEL_DATABASE_1_16_PNG_FORMAT = 1;



////////////////////////////////////////////////////////////////////////
// R2PixelDatabaseEntry class definition
////////////////////////////////////////////////////////////////////////

struct R2PixelDatabaseEntry {
public:
  R2PixelDatabaseEntry(const char *key = NULL,
    int format = 0, unsigned int size = 0, unsigned long long seek = 0,
    double offset = 0, double scale = 1, double exponent = 1)
    : format(format),
      size(size),
      seek(seek),
      offset(offset),
      scale(scale),
      exponent(exponent)
  {
    this->key[0] = '\0';
    if (key) strncpy(this->key, key, 127);
    this->key[127]='\0';
  };
  
public:
  char key[128];
  unsigned int format;
  unsigned int size;
  unsigned long long seek;
  double offset;
  double scale;
  double exponent;
};



////////////////////////////////////////////////////////////////////////
// CONSTRUCTORS/DESTRUCTORS
////////////////////////////////////////////////////////////////////////

R2PixelDatabase::
R2PixelDatabase(void)
  : fp(NULL),
    name(NULL),
    filename(NULL),
    rwaccess(NULL),
    major_version(current_major_version),
    minor_version(current_minor_version),
    swap_endian(0),
    entries_count(0),
    entries_seek(0),
    map(),
    keys()
{
}



R2PixelDatabase::
R2PixelDatabase(const R2PixelDatabase& database)
  : fp(NULL),
    name(NULL),
    filename(NULL),
    rwaccess(NULL),
    major_version(current_major_version),
    minor_version(current_minor_version),
    swap_endian(0),
    entries_count(0),
    entries_seek(0),
    map(),
    keys()
{
  RNAbort("Not implemented");
}



R2PixelDatabase::
~R2PixelDatabase(void)
{
  // Delete name
  if (name) free(name);

  // Delete filename
  if (filename) free(filename);

  // Delete rwaccess
  if (rwaccess) free(rwaccess);
}



////////////////////////////////////////////////////////////////////////
// PROPERTY MANIPULATION FUNCTIONS
////////////////////////////////////////////////////////////////////////

void R2PixelDatabase::
SetName(const char *name)
{
  // Set node name
  if (this->name) delete this->name;
  this->name = RNStrdup(name);
}
  


////////////////////////////////////////////////////////////////////////
// ENTRY ACCESS FUNCTIONS
////////////////////////////////////////////////////////////////////////

int R2PixelDatabase::
FindImage(const char *key, R2Image *image) const
{
  // Find entry
  R2PixelDatabaseEntry entry;
  if (!map.Find(key, &entry)) return FALSE;

  // Read image
  if (image) {
    // Seek to start of entry
    RNFileSeek(fp, entry.seek, RN_FILE_SEEK_SET);
  
    // Read entry
    if (!image->ReadPNGStream(fp)) {
      RNFail("Error reading %s from pixel database\n", key);
      return FALSE;
    }
  }

  // Return success
  return TRUE;
}



int R2PixelDatabase::
FindGrid(const char *key, R2Grid *grid) const
{
  // Find entry
  R2PixelDatabaseEntry entry;
  if (!map.Find(key, &entry)) return FALSE;

  // Read grid
  if (grid) {
    // Seek to start of entry
    RNFileSeek(fp, entry.seek, RN_FILE_SEEK_SET);
  
    // Read entry
    if (!grid->ReadPNGStream(fp)) {
      RNFail("Error reading %s from pixel database\n", key);
      return FALSE;
    }

    // Apply offset 
    if (entry.offset != 0) {
      grid->Subtract(entry.offset);
    }

    // Apply scale 
    if ((entry.scale != 0) && (entry.scale != 1)) {
      grid->Multiply(1.0 / entry.scale);
    }

    // Apply exponent
    if ((entry.exponent != 0) && (entry.exponent != 1)) {
      grid->Pow(1.0 / entry.exponent);
    }
  }

  // Return success
  return TRUE;
}



////////////////////////////////////////////////////////////////////////
// ENTRY MANIPULATION FUNCTIONS
////////////////////////////////////////////////////////////////////////

int R2PixelDatabase::
InsertImage(const char *key, const R2Image& image)
{
  // Seek to end of entries
  unsigned long long seek = entries_seek;
  RNFileSeek(fp, seek, RN_FILE_SEEK_SET);

  // Write pixels to file
  if (!image.WritePNGStream(fp)) return FALSE;

  // Update entries seek
  entries_seek = RNFileTell(fp);
  unsigned int size = entries_seek - seek;
  
  // Insert entry into map
  R2PixelDatabaseEntry entry(key, R2_PIXEL_DATABASE_3_8_PNG_FORMAT, size, seek);
  map.Insert(key, entry);

  // Insert key
  keys.push_back(key);
  
  // Increment number of entries
  entries_count++;

  // Return success
  return 1;
}



int R2PixelDatabase::
InsertGrid(const char *key, const R2Grid& grid,
  double offset, double scale, double exponent, RNBoolean apply_transformation)
{
  // Seek to end of entries
  unsigned long long seek = entries_seek;
  RNFileSeek(fp, seek, RN_FILE_SEEK_SET);

  // Check if need to apply offset, scale, or exponent
  if ((apply_transformation) &&
      ((offset != 0) || (scale != 1) || (exponent != 1))) {
    // Write processed pixels to file
    R2Grid tmp(grid);
    if (exponent != 1) tmp.Pow(exponent);
    if (scale != 1) tmp.Multiply(scale);
    if (offset != 0) tmp.Add(offset);
    if (!tmp.WritePNGStream(fp)) return FALSE;
  }
  else {
    // Write original pixels to file
    if (!grid.WritePNGStream(fp)) return FALSE;
  }

  // Update entries seek
  entries_seek = RNFileTell(fp);
  unsigned int size = entries_seek - seek;
  
  // Insert entry into map
  R2PixelDatabaseEntry entry(key, R2_PIXEL_DATABASE_1_16_PNG_FORMAT, size, seek, offset, scale, exponent);
  map.Insert(key, entry);
  
  // Insert key
  keys.push_back(key);
  
  // Increment number of entries
  entries_count++;

  // Return success
  return 1;
}



int R2PixelDatabase::
Remove(const char *key)
{
  // Remove entry from map
  map.Remove(key);

  // Remove key
  std::vector<std::string>::iterator it = std::find(keys.begin(), keys.end(), key);
  if (it != keys.end()) keys.erase(it);
  
  // Leaves hole in file where pixels were :(

  // Decrement number of entries
  entries_count--;

  // Return success
  return 1;
}



////////////////////////////////////////////////////////////////////////
// FILE I/O FUNCTIONS
////////////////////////////////////////////////////////////////////////

int R2PixelDatabase::
OpenFile(const char *filename, const char *rwaccess)
{
  // Remember file name
  if (this->filename) free(this->filename);
  this->filename = RNStrdup(filename);

  // Parse rwaccess
  if (this->rwaccess) free(this->rwaccess);
  if (!rwaccess) this->rwaccess = RNStrdup("w+b");
  else if (strstr(rwaccess, "w")) this->rwaccess = RNStrdup("w+b");
  else if (strstr(rwaccess, "+")) this->rwaccess = RNStrdup("r+b");
  else this->rwaccess = RNStrdup("rb"); 

  // Open file
  fp = fopen(filename, this->rwaccess);
  if (!fp) {
    RNFail("Unable to open database file %s with rwaccess %s\n", filename, rwaccess);
    return 0;
  }

  // Check if file is new
  if (!strcmp(this->rwaccess, "w+b")) {
    // Just checking ...
    assert(entries_count == 0);

    // File is new -- write header
    if (!WriteHeader(fp, 0)) return 0;

    // Update entries info
    entries_seek = RNFileTell(fp);
  }
  else {
    // Read header
    if (!ReadHeader(fp)) return 0;

    // Read entries
    if (!ReadEntries(fp, swap_endian)) return 0;
  }

  // Return success
  return 1;
}



int R2PixelDatabase::
CloseFile(void)
{
  // Check if writing file
  if (strcmp(rwaccess, "rb")) {
    // Write entries
    if (!WriteEntries(fp, swap_endian)) return 0;

    // Write header again (now that the seek values have been filled in)
    if (!WriteHeader(fp, swap_endian)) return 0;
  }

  // Close file
  fclose(fp);
  fp = NULL;
  
  // Reset filename
  if (filename) free(filename);
  filename = NULL;

  // Reset rwaccess
  if (rwaccess) free(rwaccess);
  rwaccess = NULL;

  // Return success
  return 1;
}



////////////////////////////////////////////////////////////////////////
// INTERNAL I/O FUNCTIONS
////////////////////////////////////////////////////////////////////////

int R2PixelDatabase::
ReadHeader(FILE *fp)
{
  // Seek to beginning of file
  RNFileSeek(fp, 0, RN_FILE_SEEK_SET);

  // Read magic string
  char buffer[1024]; 
  if (!RNReadChar(fp, buffer, 32, 0)) return 0;
  if (strcmp(buffer, "R2PixelDatabase")) {
    RNFail("Incorrect header (%s) in database file %s\n", buffer, filename);
    return 0;
  }

  // Read endian test
  unsigned int endian_test1, endian_test2;
  if (!RNReadUnsignedInt(fp, &endian_test1, 1, 0)) return 0;
  if (endian_test1 != 1) swap_endian = 1;
  if (!RNReadUnsignedInt(fp, &endian_test2, 1, swap_endian)) return 0;
  if (endian_test2 != 1) {
    RNFail("Incorrect endian (%x) in database file %s\n", endian_test1, filename);
    return 0;
  }

  // Read version
  if (!RNReadUnsignedInt(fp, &major_version, 1, swap_endian)) return 0;
  if (!RNReadUnsignedInt(fp, &minor_version, 1, swap_endian)) return 0;
  
  // Read entry info
  if (!RNReadUnsignedLongLong(fp, &entries_seek, 1, swap_endian)) return 0;
  if (!RNReadUnsignedInt(fp, &entries_count, 1, swap_endian)) return 0;

  // Read extra at end of header
  if (!RNReadChar(fp, buffer, 1024, swap_endian)) return 0;
  
  // Return success
  return 1;
}

 
int R2PixelDatabase::
WriteHeader(FILE *fp, int swap_endian)
{
  // Seek to beginning of file
  RNFileSeek(fp, 0, RN_FILE_SEEK_SET);

  // Get convenient variables
  unsigned int endian_test = 1;
  char magic[32] = { '\0' };
  strncpy(magic, "R2PixelDatabase", 32);
  char buffer[1024] = { '\0' };

  // Write header
  if (!RNWriteChar(fp, magic, 32, swap_endian)) return 0;
  if (!RNWriteUnsignedInt(fp, &endian_test, 1, swap_endian)) return 0;
  if (!RNWriteUnsignedInt(fp, &endian_test, 1, swap_endian)) return 0;
  if (!RNWriteUnsignedInt(fp, &major_version, 1, swap_endian)) return 0;
  if (!RNWriteUnsignedInt(fp, &minor_version, 1, swap_endian)) return 0;
  if (!RNWriteUnsignedLongLong(fp, &entries_seek, 1, swap_endian)) return 0;
  if (!RNWriteUnsignedInt(fp, &entries_count, 1, swap_endian)) return 0;
  if (!RNWriteChar(fp, buffer, 1024, swap_endian)) return 0;
  
  // Return success
  return 1;
}


 
int R2PixelDatabase::
ReadEntries(FILE *fp, int swap_endian)
{
  // Seek to entries seek
  RNFileSeek(fp, entries_seek, RN_FILE_SEEK_SET);

  // Read entries
  int dummy = 0;
  for (unsigned int i = 0; i < entries_count; i++) {
    R2PixelDatabaseEntry entry;
    if (!RNReadChar(fp, entry.key, 128, swap_endian)) return 0;
    if (!RNReadUnsignedInt(fp, &entry.format, 1, swap_endian)) return 0;
    if (!RNReadUnsignedInt(fp, &entry.size, 1, swap_endian)) return 0;
    if (!RNReadUnsignedLongLong(fp, &entry.seek, 1, swap_endian)) return 0;
    if (!RNReadDouble(fp, &entry.offset, 1, swap_endian)) return 0;
    if (!RNReadDouble(fp, &entry.scale, 1, swap_endian)) return 0;
    if (!RNReadDouble(fp, &entry.exponent, 1, swap_endian)) return 0;
    for (int j = 0; j < 12; j++) RNReadInt(fp, &dummy, 1, swap_endian);
    map.Insert(entry.key, entry);
    keys.push_back(entry.key);
  }

  // Return success
  return 1;
}


 
int R2PixelDatabase::
WriteEntries(FILE *fp, int swap_endian)
{
  // Seek to entries seek
  RNFileSeek(fp, entries_seek, RN_FILE_SEEK_SET);

  // Write entries
  int dummy = 0;
  char buffer[128] = { '\0' };
  std::map<std::string, R2PixelDatabaseEntry, RNMapComparator<std::string> >::iterator it;
  for (it = map.m->begin(); it != map.m->end(); ++it) {
    R2PixelDatabaseEntry entry = it->second;
    strncpy(buffer, entry.key, 128);
    if (!RNWriteChar(fp, buffer, 128, swap_endian)) return 0;
    if (!RNWriteUnsignedInt(fp, &entry.format, 1, swap_endian)) return 0;
    if (!RNWriteUnsignedInt(fp, &entry.size, 1, swap_endian)) return 0;
    if (!RNWriteUnsignedLongLong(fp, &entry.seek, 1, swap_endian)) return 0;
    if (!RNWriteDouble(fp, &entry.offset, 1, swap_endian)) return 0;
    if (!RNWriteDouble(fp, &entry.scale, 1, swap_endian)) return 0;
    if (!RNWriteDouble(fp, &entry.exponent, 1, swap_endian)) return 0;
    for (int j = 0; j < 12; j++) RNWriteInt(fp, &dummy, 1, swap_endian);
  }

  // Return success
  return 1;
}

 
 
} // namespace gaps
