/* Source file for the R3 surfel database class */



////////////////////////////////////////////////////////////////////////
// INCLUDE FILES
////////////////////////////////////////////////////////////////////////

#include "R3Surfels.h"



////////////////////////////////////////////////////////////////////////
// Namespace
////////////////////////////////////////////////////////////////////////

namespace gaps {



////////////////////////////////////////////////////////////////////////
// PRINT DEBUG CONTROL
////////////////////////////////////////////////////////////////////////

// #define PRINT_DEBUG



////////////////////////////////////////////////////////////////////////
// Versioning variables
////////////////////////////////////////////////////////////////////////

static unsigned int current_major_version = 4;
static unsigned int current_minor_version = 1;



////////////////////////////////////////////////////////////////////////
// CONSTRUCTORS/DESTRUCTORS
////////////////////////////////////////////////////////////////////////

R3SurfelDatabase::
R3SurfelDatabase(void)
  : fp(NULL),
    filename(NULL),
    rwaccess(NULL),
    major_version(current_major_version),
    minor_version(current_minor_version),
    swap_endian(0),
    file_blocks_offset(0),
    file_blocks_count(0),
    blocks(),
    nsurfels(0),
    bbox(FLT_MAX,FLT_MAX,FLT_MAX,-FLT_MAX,-FLT_MAX,-FLT_MAX),
    timestamp_range(FLT_MAX,-FLT_MAX),
    max_identifier(0),
    name(NULL),
    tree(NULL),
    resident_surfels(0)
{
}



R3SurfelDatabase::
R3SurfelDatabase(const R3SurfelDatabase& database)
  : fp(NULL),
    filename(NULL),
    rwaccess(NULL),
    major_version(current_major_version),
    minor_version(current_minor_version),
    swap_endian(0),
    file_blocks_offset(0),
    file_blocks_count(0),
    blocks(),
    nsurfels(0),
    bbox(FLT_MAX,FLT_MAX,FLT_MAX,-FLT_MAX,-FLT_MAX,-FLT_MAX),
    timestamp_range(FLT_MAX,-FLT_MAX),
    max_identifier(0),
    name(RNStrdup(database.name)),
    tree(NULL),
    resident_surfels(0)
{
  RNAbort("Not implemented");
}



R3SurfelDatabase::
~R3SurfelDatabase(void)
{
  // Close database
  if (IsOpen()) CloseFile();

  // Delete blocks
  while (NBlocks() > 0) {
    R3SurfelBlock *block = Block(NBlocks()-1);
    // This allows delete of database without releasing all blocks
    block->file_read_count = 0; 
    RemoveBlock(block);
    delete block;
  }

  // Remove from tree
  if (tree) tree->database = NULL;

  // Delete filename
  if (filename) free(filename);

  // Delete rwaccess
  if (rwaccess) free(rwaccess);

  // Delete name
  if (name) free(name);
}



////////////////////////////////////////////////////////////////////////
// PROPERTY MANIPULATION FUNCTIONS
////////////////////////////////////////////////////////////////////////

void R3SurfelDatabase::
SetName(const char *name)
{
  // Set node name
  if (this->name) delete this->name;
  this->name = RNStrdup(name);
}



////////////////////////////////////////////////////////////////////////
// SURFEL MANIPULATION FUNCTIONS
////////////////////////////////////////////////////////////////////////

void R3SurfelDatabase::
SetMarks(RNBoolean mark)
{
  // Mark all surfels
  for (int i = 0; i < NBlocks(); i++) {
    R3SurfelBlock *block = Block(i);
    block->SetMarks(mark);
  }
}



////////////////////////////////////////////////////////////////////////
// BLOCK MANIPULATION FUNCTIONS
////////////////////////////////////////////////////////////////////////

void R3SurfelDatabase::
InsertBlock(R3SurfelBlock *block)
{
  // Just checking
  assert(block->database == NULL);
  assert(block->database_index == -1);
  assert(block->file_surfels_offset == 0);
  assert(block->file_surfels_count == 0);
  assert(block->file_read_count == 0);

  // Update block database info
  block->database = this;
  block->database_index = blocks.NEntries();
  block->file_surfels_offset = 0;
  block->file_surfels_count = 0;
  block->file_read_count = (block->surfels) ? 1 : 0;
  block->SetDirty(TRUE);

  // Insert block
  blocks.Insert(block);

  // Update bounding box
  bbox.Union(block->BBox());

  // Update timestamp range
  timestamp_range.Union(block->TimestampRange());

  // Update max identifier
  if (block->MaxIdentifier() > max_identifier)
    max_identifier = block->MaxIdentifier();

  // Update number of surfels
  nsurfels += block->NSurfels();

  // Update block
  block->UpdateAfterInsert(this);

  // Update resident surfels
  if (block->surfels) resident_surfels += block->NSurfels();

#ifdef PRINT_DEBUG
  // Print debug message
  printf("Inserted Block %6d : %6d %9ld : %9.3f %9.3f %9.3f\n", 
    block->database_index, block->nsurfels, resident_surfels,
    block->Centroid().X(), block->Centroid().Y(), block->Centroid().Z()); 
  fflush(stdout);
#endif
}



void R3SurfelDatabase::
RemoveBlock(R3SurfelBlock *block)
{
  // Just checking
  assert(block->file_read_count == 0);
  assert(block->database == this);
  assert(block->node == NULL);
    
  // Update resident surfels
  if (block->surfels) resident_surfels -= block->NSurfels();
  assert(resident_surfels >= 0);
    
  // Update block
  block->UpdateBeforeRemove(this);
    
  // Find block
  RNArrayEntry *entry = blocks.KthEntry(block->database_index);
  if (!entry) return;
  R3SurfelBlock *tail = blocks.Tail();
  blocks.EntryContents(entry) = tail;
  tail->database_index = block->database_index;
  blocks.RemoveTail();
    
  // Reset block database info
  block->database = NULL;
  block->database_index = -1;
  block->file_surfels_offset = 0;
  block->file_surfels_count = 0;
  block->file_read_count = 0;
  block->SetDirty(FALSE);
    
  // Update number of surfels
  nsurfels -= block->NSurfels();
  assert(nsurfels >= 0);

  // Does not update bounding box
  // XXX

#ifdef PRINT_DEBUG
  // Print debug message
  printf("Removed Block  %6d : %6d %9ld : %9.3f %9.3f %9.3f\n", 
         block->database_index, block->nsurfels, resident_surfels,
         block->Centroid().X(), block->Centroid().Y(), block->Centroid().Z()); 
  fflush(stdout);
#endif
}



void R3SurfelDatabase::
RemoveAndDeleteBlock(R3SurfelBlock *block)
{
  // Check if still referenced
  if (block->file_read_count == 0) {
    // Block is not referenced, can simply delete it
    RemoveBlock(block);
    delete block;
  }
  else {
    // Block is referenced, mark for delete later
    block->flags.Add(R3_SURFEL_BLOCK_DELETE_PENDING_FLAG);
  }
}



////////////////////////////////////////////////////////////////////////
// HIGH-LEVEL MANIPULATION FUNCTIONS
////////////////////////////////////////////////////////////////////////

int R3SurfelDatabase::
InsertSubsetBlocks(R3SurfelBlock *block, 
  const RNArray<const R3Surfel *>& subset1, const RNArray<const R3Surfel *>& subset2, 
  R3SurfelBlock **blockA, R3SurfelBlock **blockB)
{
  // Check subset1
  if (subset1.IsEmpty()) {
    if (blockA) *blockA = NULL;
    if (blockB) *blockB = block;
    return 0;
  }

  // Check subset2
  if (subset2.IsEmpty()) {
    if (blockA) *blockA = block;
    if (blockB) *blockB = NULL;
    return 0;
  }

  // Create new blocks
  R3SurfelBlock *block1 = new R3SurfelBlock(subset1, block->PositionOrigin(), block->TimestampOrigin());
  R3SurfelBlock *block2 = new R3SurfelBlock(subset2, block->PositionOrigin(), block->TimestampOrigin());
    
  // Insert new blocks 
  InsertBlock(block1);
  InsertBlock(block2);

  // Update file offsets
  if ((block->file_surfels_offset > 0) && (block->file_surfels_count > 0)) {
    block1->file_surfels_offset = block->file_surfels_offset;
    block1->file_surfels_count = block1->NSurfels();
    block2->file_surfels_offset = block->file_surfels_offset + block1->NSurfels() * sizeof(R3Surfel);
    block2->file_surfels_count = block2->NSurfels();
    block->file_surfels_offset = 0;
    block->file_surfels_count = 0;
  }

  // Update file read counts ???
  if (block->file_read_count > 0) {
    block1->file_read_count = block->file_read_count;
    block2->file_read_count = block->file_read_count;
  }
    
  // Update block properties
  block1->UpdateProperties();
  block2->UpdateProperties();

  // Release blocks
  ReleaseBlock(block1);
  ReleaseBlock(block2);
      
  // Return new blocks
  if (blockA) *blockA = block1;
  if (blockB) *blockB = block2;

  // Return success
  return 1;
}
  
  

////////////////////////////////////////////////////////////////////////
// DISPLAY FUNCTIONS
////////////////////////////////////////////////////////////////////////

void R3SurfelDatabase::
Draw(RNFlags flags) const
{
  // Draw blocks
  for (int i = 0; i < blocks.NEntries(); i++) {
    blocks[i]->Draw(flags);
  }
}



void R3SurfelDatabase::
Print(FILE *fp, const char *prefix, const char *suffix) const
{
  // Check fp
  if (!fp) fp = stdout;

  // Print name
  fprintf(fp, "%s%s%s\n", (prefix) ? prefix : "", (name) ? name : "Database", (suffix) ? suffix : "");

  // Add indent to prefix
  char indent_prefix[1024];
  sprintf(indent_prefix, "%s  ", (prefix) ? prefix : "");

  // Print blocks
  for (int i = 0; i < NBlocks(); i++) {
    blocks[i]->Print(fp, indent_prefix, suffix);
  }
}



////////////////////////////////////////////////////////////////////////
// I/O FUNCTIONS
////////////////////////////////////////////////////////////////////////

int R3SurfelDatabase::
ReadFile(const char *filename)
{
  // Parse input filename extension
  const char *extension;
  if (!(extension = strrchr(filename, '.'))) {
    printf("Filename %s has no extension (e.g., .ply)\n", filename);
    return 0;
  }

  // Read file of appropriate type
  if (!strncmp(extension, ".list", 5)) {
    // Open file
    FILE *fp;
    if (!(fp = fopen(filename, "r"))) {
      RNFail("Unable to open file %s\n", filename);
      return 0;
    }

    // Read file names
    char buffer[4096];
    while (fscanf(fp, "%s", buffer) == (unsigned int) 1) {
      // Create block
      R3SurfelBlock *block = new R3SurfelBlock();
      if (!block) {
        RNFail("Unable to create block for %s\n", buffer);
        return 0;
      }
    
      // Read block
      if (!block->ReadFile(buffer)) { 
        delete block; 
        return 0; 
      }

      // Update properties
      block->UpdateProperties();

      // Insert block
      InsertBlock(block);

      // Release block
      ReleaseBlock(block);
    }

    // Close file
    fclose(fp);
  }
  else { 
    // Create block
    R3SurfelBlock *block = new R3SurfelBlock();
    if (!block) {
      RNFail("Unable to create block\n");
      return 0;
    }
    
    // Read file
    if (!block->ReadFile(filename)) { 
      delete block; 
      return 0; 
    }
    
    // Update properties
    block->UpdateProperties();

    // Insert block
    InsertBlock(block);

    // Return success
    return 1;
  }

  // Should never get here
  return 0;
}



int R3SurfelDatabase::
WriteFile(const char *filename) const
{
  // Parse input filename extension
  const char *extension;
  if (!(extension = strrchr(filename, '.'))) {
    printf("Filename %s has no extension (e.g., .xyz)\n", filename);
    return 0;
  }

  // Write file of appropriate type
  if (NBlocks() == 1) {
    R3SurfelDatabase *database = (R3SurfelDatabase *) this;
    R3SurfelBlock *block = Block(0);
    if (!database->ReadBlock(block)) return 0;
    if (!block->WriteFile(filename)) return 0;
    if (!database->ReleaseBlock(block)) return 0;
  }
  else {
    RNFail("Invalid file extension %s for database with more than one block\n", extension);
    return 0;
  }

  // Return success
  return 1;
}



////////////////////////////////////////////////////////////////////////
// I/O UTILITY FUNCTIONS
////////////////////////////////////////////////////////////////////////

int R3SurfelDatabase::
ReadSurfel(FILE *fp, R3Surfel *ptr, int count, int swap_endian,
  unsigned int major_version, unsigned int minor_version) const
{
  // Check database version
  if ((major_version == current_major_version) && (minor_version == current_minor_version)) {
    int sofar = 0;
    while (sofar < count) {
      size_t status = fread(ptr, sizeof(R3Surfel), count - sofar, fp);
      if (status > 0) sofar += status;
      else { RNFail("Unable to read surfel from database file\n"); return 0; }
    }
  }
  else {
    if (major_version == 3) {
      for (int i = 0; i < count; i++) {
        fread(ptr[i].position, sizeof(float), 3, fp);
        fread(ptr[i].normal, sizeof(RNInt16), 3, fp);
        fread(ptr[i].radius, sizeof(RNInt16), 1, fp);
        fread(ptr[i].color, sizeof(unsigned char), 3, fp);
        fread(&ptr[i].flags, sizeof(unsigned char), 1, fp);
      }
    }
    else if (major_version < 2) {
      for (int i = 0; i < count; i++) {
        float position[3];
        unsigned char color_and_flags[4];
        fread(position, sizeof(float), 3, fp);
        fread(color_and_flags, sizeof(unsigned char), 4, fp);
        ptr[i].SetPosition(position);
        ptr[i].SetColor(color_and_flags);
        ptr[i].SetFlags(color_and_flags[3]);
      }
    }
  }

  // Swap endian
  if (swap_endian) {
    for (int i = 0; i < count; i++) {
      RNSwap4(ptr[i].position, 3);
      RNSwap2(ptr[i].normal, 3);
      RNSwap2(ptr[i].tangent, 3);
      RNSwap2(ptr[i].radius, 2);
      RNSwap4(&ptr[i].identifier, 1);
    }
  }

  // Return success
  return 1;
}



int R3SurfelDatabase::
WriteSurfel(FILE *fp, R3Surfel *ptr, int count, int swap_endian, 
  unsigned int major_version, unsigned int minor_version) const
{
  // Swap endian
  if (swap_endian) {
    for (int i = 0; i < count; i++) {
      RNSwap4(ptr[i].position, 3);
      RNSwap2(ptr[i].normal, 3);
      RNSwap2(ptr[i].tangent, 3);
      RNSwap2(ptr[i].radius, 2);
      RNSwap4(&ptr[i].identifier, 1);
    }
  }

  // Clear surfel marks
  for (int i = 0; i < count; i++) ptr[i].SetMark(FALSE);

  // Write current version of surfel
  int sofar = 0;
  int status = 1;
  while (sofar < count) {
    size_t n = fwrite(ptr, sizeof(R3Surfel), count - sofar, fp);
    if (n > 0) sofar += n;
    else { RNFail("Unable to write surfel to database file\n"); status = 0; }
  }

  // Swap endian back
  if (swap_endian) {
    for (int i = 0; i < count; i++) {
      RNSwap4(ptr[i].position, 3);
      RNSwap2(ptr[i].normal, 3);
      RNSwap2(ptr[i].tangent, 3);
      RNSwap2(ptr[i].radius, 2);
      RNSwap4(&ptr[i].identifier, 1);
    }
  }

  // Return status
  return status;
}



////////////////////////////////////////////////////////////////////////
// BLOCK I/O FUNCTIONS
////////////////////////////////////////////////////////////////////////

int R3SurfelDatabase::
InternalReadBlock(R3SurfelBlock *block)
{
  // Check number of surfels
  if (block->NSurfels() == 0) return 1;

  // Just checking
  assert(fp);
  assert(block->database == this);
  assert(block->file_surfels_offset > 0);
  assert(block->file_surfels_count >= (unsigned int) block->nsurfels);

  // Allocate surfels
  block->surfels = new R3Surfel [ block->nsurfels ];
  if (!block->surfels) {
    RNFail("Unable to allocate surfels\n");
    return 0;
  }
  
  // Read surfels
  RNFileSeek(fp, block->file_surfels_offset, RN_FILE_SEEK_SET);
  if (!ReadSurfel(fp, block->surfels, block->nsurfels, swap_endian, major_version, minor_version)) return 0;
  
  // Update resident surfels
  resident_surfels += block->NSurfels();

#ifdef PRINT_DEBUG
  // Print debug message
  printf("Read Block     %6d : %6d %9ld : %9.3f %9.3f %9.3f\n", 
         block->database_index, block->nsurfels, resident_surfels,
         block->Centroid().X(), block->Centroid().Y(), block->Centroid().Z()); 
  fflush(stdout);
#endif

  // Return success
  return 1;
}



int R3SurfelDatabase::
InternalReleaseBlock(R3SurfelBlock *block)
{
  // Just checking
  assert(block->database == this);

  // Write block
  if (!SyncBlock(block)) return 0;
    
#ifdef R3_SURFEL_DRAW_WITH_DISPLAY_LIST
  // Delete opengl display lists
  if (block->opengl_id > 0) {
    glDeleteLists(block->opengl_id, 2);
    block->opengl_id = 0;
  }
#endif

#ifdef DRAW_WITH_VBO
  // Delete opengl vertex buffer object
  if (block->opengl_id > 0) {
    glDeleteBuffers(1, &block->opengl_id);
    block->opengl_id = 0;
  }
#endif

  // Delete surfels
  if (block->surfels) {
    delete [] block->surfels;
    block->surfels = NULL;
  }
      
  // Update resident surfels
  resident_surfels -= block->NSurfels();
  assert(resident_surfels >= 0);

#ifdef PRINT_DEBUG
  // Print debug message
  printf("Released Block %6d : %6d %9ld : %9.3f %9.3f %9.3f\n", 
         block->database_index, block->nsurfels, resident_surfels,
         block->Centroid().X(), block->Centroid().Y(), block->Centroid().Z()); 
  fflush(stdout);
#endif

  // Return success
  return 1;
}



int R3SurfelDatabase::
InternalSyncBlock(R3SurfelBlock *block)
{
  // Check file
  if (!fp) return 1;

  // Check number of surfels
  if (block->NSurfels() == 0) return 1;

  // Check file rwaccess
  if (!strstr(rwaccess, "+")) {
    RNFail("Unable to write block to read-only file\n");
    return 0;
  }

  // Check database version
  if ((major_version != current_major_version) || (minor_version != current_minor_version)) {
    RNFail("Unable to write block to database with different version\n");
    return 0;
  }

  // Just checking
  assert(block->database == this);

  // Check if surfels can be put at original offset in file
  if ((block->file_surfels_offset > 0) && ((unsigned int) block->nsurfels <= block->file_surfels_count)) {
    // Surfels fit at original offset in file
    RNFileSeek(fp, block->file_surfels_offset, RN_FILE_SEEK_SET);
  }
  else {
    // Surfels must be put at end of file
    RNFileSeek(fp, 0, RN_FILE_SEEK_END);
    block->file_surfels_offset = RNFileTell(fp);
    block->file_surfels_count = block->nsurfels;
  }

  // Write surfels to file
  if (!WriteSurfel(fp, block->surfels, block->nsurfels, swap_endian, major_version, minor_version)) return 0;

#ifdef PRINT_DEBUG
  // Print debug message
  printf("Synced Block %6d : %6d %9ld : %9.3f %9.3f %9.3f\n", 
         block->database_index, block->nsurfels, resident_surfels,
         block->Centroid().X(), block->Centroid().Y(), block->Centroid().Z()); 
  fflush(stdout);
#endif

  // Return success
  return 1;
}



////////////////////////////////////////////////////////////////////////
// FILE I/O FUNCTIONS
////////////////////////////////////////////////////////////////////////

int R3SurfelDatabase::
WriteHeader(FILE *fp, int swap_endian)
{
  // Get convenient variables
  unsigned int endian_test = 1;
  unsigned int nblocks = blocks.NEntries();
  char magic[32] = { '\0' };
  strncpy(magic, "R3SurfelDatabase", 32);
  char buffer[1024] = { '\0' };

  // Write header
  RNFileSeek(fp, 0, RN_FILE_SEEK_SET);
  if (!RNWriteChar(fp, magic, 32, swap_endian)) return 0;
  if (!RNWriteUnsignedInt(fp, &endian_test, 1, swap_endian)) return 0;
  if (!RNWriteUnsignedInt(fp, &endian_test, 1, swap_endian)) return 0;
  if (!RNWriteUnsignedInt(fp, &major_version, 1, swap_endian)) return 0;
  if (!RNWriteUnsignedInt(fp, &minor_version, 1, swap_endian)) return 0;
  if (!RNWriteUnsignedLongLong(fp, &file_blocks_offset, 1, swap_endian)) return 0;
  if (!RNWriteUnsignedInt(fp, &file_blocks_count, 1, swap_endian)) return 0;
  if (!RNWriteUnsignedInt(fp, &nblocks, 1, swap_endian)) return 0;
  if (!RNWriteLongLong(fp, &nsurfels, 1, swap_endian)) return 0;
  if (!RNWriteDouble(fp, &bbox[0][0], 6, swap_endian)) return 0;
  if (!RNWriteDouble(fp, &timestamp_range[0], 2, swap_endian)) return 0;
  if (!RNWriteUnsignedInt(fp, &max_identifier, 1, swap_endian)) return 0;
  if (!RNWriteChar(fp, buffer, 1004, swap_endian)) return 0;
  
  // Return success
  return 1;
}



int R3SurfelDatabase::
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
    // File is new -- write header
    if (!WriteHeader(fp, 0)) return 0;
  }
  else {
    // Read unique string
    char buffer[1024]; 
    if (!RNReadChar(fp, buffer, 32, 0)) return 0;
    if (strcmp(buffer, "R3SurfelDatabase")) {
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
    if ((major_version < 2) || (major_version > 4)) {
      RNFail("Incorrect version (%d.%d) in database file %s\n", major_version, minor_version, filename);
      return 0;
    }
  
    // Read block info
    unsigned int nblocks;
    if (!RNReadUnsignedLongLong(fp, &file_blocks_offset, 1, swap_endian)) return 0;
    if (!RNReadUnsignedInt(fp, &file_blocks_count, 1, swap_endian)) return 0;
    if (!RNReadUnsignedInt(fp, &nblocks, 1, swap_endian)) return 0;

    // Read number of surfels
    if (major_version < 4) {
      // nsurfels used to be an int
      int nsurfels32;
      if (!RNReadInt(fp, &nsurfels32, 1, swap_endian)) return 0;
      nsurfels = nsurfels32;
    }
    else {
      // nsurfels now is an RNInt64
      if (!RNReadLongLong(fp, &nsurfels, 1, swap_endian)) return 0;
    }

    // Read bounding box
    if (!RNReadDouble(fp, &bbox[0][0], 6, swap_endian)) return 0;

    // Read timestamp range
    if (!RNReadDouble(fp, &timestamp_range[0], 2, swap_endian)) return 0;

    // Read max identifier
    if (!RNReadUnsignedInt(fp, &max_identifier, 1, swap_endian)) return 0;

    // Read extra at end of header
    if (!RNReadChar(fp, buffer, 1004, swap_endian)) return 0;

    // Read blocks
    RNFileSeek(fp, file_blocks_offset, RN_FILE_SEEK_SET);
    for (unsigned int i = 0; i < nblocks; i++) {
      R3SurfelBlock *block = new R3SurfelBlock();
      unsigned int block_flags;
      if (!RNReadUnsignedLongLong(fp, &block->file_surfels_offset, 1, swap_endian)) return 0;
      if (!RNReadUnsignedInt(fp, &block->file_surfels_count, 1, swap_endian)) return 0;
      if (!RNReadInt(fp, &block->nsurfels, 1, swap_endian)) return 0;
      if (!RNReadDouble(fp, &block->position_origin[0], 3, swap_endian)) return 0;
      if (!RNReadDouble(fp, &block->bbox[0][0], 6, swap_endian)) return 0;
      if (!RNReadDouble(fp, &block->resolution, 1, swap_endian)) return 0;
      if (!RNReadUnsignedInt(fp, &block_flags, 1, swap_endian)) return 0;
      if (!RNReadDouble(fp, &block->timestamp_origin, 1, swap_endian)) return 0;
      if (!RNReadDouble(fp, &block->timestamp_range[0], 2, swap_endian)) return 0;
      if (!RNReadUnsignedInt(fp, &block->max_identifier, 1, swap_endian)) return 0;
      if (!RNReadUnsignedInt(fp, &block->min_identifier, 1, swap_endian)) return 0;
      if (!RNReadChar(fp, buffer, 32, swap_endian)) return 0;
      block->flags = block_flags;
      block->SetDirty(FALSE);
      block->database = this;
      block->database_index = blocks.NEntries();
      blocks.Insert(block);
    }
  }

  // Return success
  return 1;
}



int R3SurfelDatabase::
SyncFile(void)
{
  // Check if file is read-only
  if (!strcmp(rwaccess, "rb")) return 1;

  // Sync blocks
  for (int i = 0; i < blocks.NEntries(); i++) {
    R3SurfelBlock *block = blocks.Kth(i);
    if (!SyncBlock(block)) return 0;
  }

  // Update blocks offset
  unsigned int nblocks = blocks.NEntries();
  if (nblocks  > file_blocks_count) {
    RNFileSeek(fp, 0, RN_FILE_SEEK_END);
    file_blocks_offset = RNFileTell(fp);
    file_blocks_count = nblocks;
  }

  // Write blocks
  char buffer[128] = { '\0' };
  RNFileSeek(fp, file_blocks_offset, RN_FILE_SEEK_SET);
  for (int i = 0; i < blocks.NEntries(); i++) {
    R3SurfelBlock *block = blocks.Kth(i);
    unsigned int block_flags = block->flags; 
    if (!RNWriteUnsignedLongLong(fp, &block->file_surfels_offset, 1, swap_endian)) return 0;
    if (!RNWriteUnsignedInt(fp, &block->file_surfels_count, 1, swap_endian)) return 0;
    if (!RNWriteInt(fp, &block->nsurfels, 1, swap_endian)) return 0;
    if (!RNWriteDouble(fp, &block->position_origin[0], 3, swap_endian)) return 0;
    if (!RNWriteDouble(fp, &block->bbox[0][0], 6, swap_endian)) return 0;
    if (!RNWriteDouble(fp, &block->resolution, 1, swap_endian)) return 0;
    if (!RNWriteUnsignedInt(fp, &block_flags, 1, swap_endian)) return 0;
    if (!RNWriteDouble(fp, &block->timestamp_origin, 1, swap_endian)) return 0;
    if (!RNWriteDouble(fp, &block->timestamp_range[0], 2, swap_endian)) return 0;
    if (!RNWriteUnsignedInt(fp, &block->max_identifier, 1, swap_endian)) return 0;
    if (!RNWriteUnsignedInt(fp, &block->min_identifier, 1, swap_endian)) return 0;
    if (!RNWriteChar(fp, buffer, 32, swap_endian)) return 0;
  }

  // Write header again (now that the offset values have been filled in)
  if (!WriteHeader(fp, swap_endian)) return 0;

  // Return success
  return 1;
}



int R3SurfelDatabase::
CloseFile(void)
{
  // Sync file
  if (!SyncFile()) return 0;

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



} // namespace gaps
