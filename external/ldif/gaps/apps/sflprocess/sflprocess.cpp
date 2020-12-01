// Source file for the surfel scene processing program



////////////////////////////////////////////////////////////////////////
// Include files
////////////////////////////////////////////////////////////////////////

namespace gaps {}
using namespace gaps;
#include "R3Surfels/R3Surfels.h"



////////////////////////////////////////////////////////////////////////
// Program arguments
////////////////////////////////////////////////////////////////////////

static char *scene_name = NULL;
static char *database_name = NULL;
static int aerial_only = 0;
static int terrestrial_only = 0;
static int print_verbose = 0;
static int print_debug = 0;



////////////////////////////////////////////////////////////////////////
// Scene I/O Functions
////////////////////////////////////////////////////////////////////////

static R3SurfelScene *
OpenScene(const char *scene_name, const char *database_name)
{
  // Start statistics
  RNTime start_time;
  start_time.Read();

  // Allocate scene
  R3SurfelScene *scene = new R3SurfelScene();
  if (!scene) {
    RNFail("Unable to allocate scene\n");
    return NULL;
  }

  // Open scene files
  if (!scene->OpenFile(scene_name, database_name, "r+", "r+")) {
    delete scene;
    return NULL;
  }

  // Print statistics
  if (print_verbose) {
    printf("Opened scene ...\n");
    printf("  Time = %.2f seconds\n", start_time.Elapsed());
    printf("  # Objects = %d\n", scene->NObjects());
    printf("  # Labels = %d\n", scene->NLabels());
    printf("  # Assignments = %d\n", scene->NLabelAssignments());
    printf("  # Features = %d\n", scene->NFeatures());
    printf("  # Scans = %d\n", scene->NScans());
    printf("  # Images = %d\n", scene->NImages());
    printf("  # Nodes = %d\n", scene->Tree()->NNodes());
    printf("  # Blocks = %d\n", scene->Tree()->Database()->NBlocks());
    printf("  # Surfels = %lld\n", scene->Tree()->Database()->NSurfels());
    fflush(stdout);
  }

  // Return scene
  return scene;
}



static int
CloseScene(R3SurfelScene *scene)
{
  // Start statistics
  RNTime start_time;
  start_time.Read();

  // Close scene files
  if (!scene->CloseFile()) return 0;

  // Print statistics
  if (print_verbose) {
    printf("Closed scene ...\n");
    printf("  Time = %.2f seconds\n", start_time.Elapsed());
    printf("  # Objects = %d\n", scene->NObjects());
    printf("  # Labels = %d\n", scene->NLabels());
    printf("  # Assignments = %d\n", scene->NLabelAssignments());
    printf("  # Features = %d\n", scene->NFeatures());
    printf("  # Scans = %d\n", scene->NScans());
    printf("  # Images = %d\n", scene->NImages());
    printf("  # Nodes = %d\n", scene->Tree()->NNodes());
    printf("  # Blocks = %d\n", scene->Tree()->Database()->NBlocks());
    printf("  # Surfels = %lld\n", scene->Tree()->Database()->NSurfels());
    fflush(stdout);
  }

  // Return success
  return 1;
}



////////////////////////////////////////////////////////////////////////
// Other I/O Functions
////////////////////////////////////////////////////////////////////////

static R2Grid *
ReadGrid(const char *filename)
{
  // Start statistics
  RNTime start_time;
  start_time.Read();

  // Allocate a grid
  R2Grid *grid = new R2Grid();
  if (!grid) {
    RNFail("Unable to allocate grid");
    return NULL;
  }

  // Read grid
  int status = grid->Read(filename);
  if (!status) {
    RNFail("Unable to read grid file %s", filename);
    return NULL;
  }

  // Print statistics
  if (print_verbose) {
    printf("Read grid from %s\n", filename);
    printf("  Time = %.2f seconds\n", start_time.Elapsed());
    printf("  Resolution = %d %d\n", grid->XResolution(), grid->YResolution());
    printf("  Spacing = %g\n", grid->GridToWorldScaleFactor());
    printf("  Cardinality = %d\n", grid->Cardinality());
    RNInterval grid_range = grid->Range();
    printf("  Minimum = %g\n", grid_range.Min());
    printf("  Maximum = %g\n", grid_range.Max());
    printf("  L1Norm = %g\n", grid->L1Norm());
    printf("  L2Norm = %g\n", grid->L2Norm());
    fflush(stdout);
  }

  // Return success
  return grid;
}



////////////////////////////////////////////////////////////////////////
// CREATE FUNCTIONS (creates structure without loading surfels)
////////////////////////////////////////////////////////////////////////

static R3SurfelNode *
CreateNode(R3SurfelScene *scene, const char *node_name, const char *parent_name)
{
  // Get surfel tree
  R3SurfelTree *tree = scene->Tree();
  if (!tree) {
    RNFail("Scene has no tree\n");
    return NULL;
  }    

  // Find parent node
  R3SurfelNode *parent_node = tree->FindNodeByName(parent_name);
  if (!parent_node) {
    RNFail("Unable to find parent node with name %s\n", parent_name);
    return NULL;
  }

  // Create node
  R3SurfelNode *node = new R3SurfelNode(node_name);
  if (!node) {
    RNFail("Unable to allocate node\n");
    return NULL;
  }
            
  // Insert node into tree
  tree->InsertNode(node, parent_node);

  // Return node
  return node;
}



static R3SurfelObject *
CreateObject(R3SurfelScene *scene, const char *object_name, 
  const char *parent_name, const char *node_name)
{
  // Get surfel tree
  R3SurfelTree *tree = scene->Tree();
  if (!tree) {
    RNFail("Scene has no tree\n");
    return NULL;
  }    

  // Find parent object
  R3SurfelObject *parent_object = scene->FindObjectByName(parent_name);
  if (!parent_object) {
    RNFail("Unable to find parent object with name %s\n", parent_name);
    return NULL;
  }

  // Find node
  R3SurfelNode *node = NULL;
  if (strcmp(node_name, "None") && strcmp(node_name, "none") && strcmp(node_name, "NONE")) {
    node = tree->FindNodeByName(node_name);
    if (!node) {
      RNFail("Unable to find parent node with name %s\n", node_name);
      return NULL;
    }
  }

  // Create object
  R3SurfelObject *object = new R3SurfelObject(object_name);
  if (!object) {
    RNFail("Unable to allocate object\n");
    return NULL;
  }
         
  // Insert node into object
  if (node) object->InsertNode(node);
  
  // Insert object into scene
  scene->InsertObject(object, parent_object);

  // Return object
  return object;
}



static R3SurfelLabel *
CreateLabel(R3SurfelScene *scene, const char *label_name, const char *parent_name)
{
  // Find parent label
  R3SurfelLabel *parent_label = scene->FindLabelByName(parent_name);
  if (!parent_label) {
    RNFail("Unable to find parent label with name %s\n", parent_name);
    return NULL;
  }

  // Create label
  R3SurfelLabel *label = CreateLabel(scene, parent_label, label_name);
  if (!label) {
    RNFail("Unable to create label\n");
    return NULL;
  }

  // Return label
  return label;
}



////////////////////////////////////////////////////////////////////////
// LOAD FUNCTIONS
////////////////////////////////////////////////////////////////////////

#if 0
static R3SurfelNode *
LoadSurfels(R3SurfelScene *scene, R3SurfelNode *parent_node, 
  R3Point *points, int npoints, const char *node_name)
{
  // Get surfel tree
  R3SurfelTree *tree = scene->Tree();
  if (!tree) {
    RNFail("Scene has no tree\n");
    return NULL;
  }    

  // Get surfel database
  R3SurfelDatabase *database = tree->Database();
  if (!database) {
    RNFail("Scene has no database\n");
    return NULL;
  }    

  // Create node
  R3SurfelNode *node = new R3SurfelNode(node_name);
  if (!node) {
    RNFail("Unable to allocate node\n");
    return NULL;
  }
            
  // Insert node into tree
  tree->InsertNode(node, parent_node);
          
  // Create block
  R3SurfelBlock *block = new R3SurfelBlock(points, npoints);
  if (!block) {
    RNFail("Unable to allocate block\n");
    return NULL;
  }
          
  // Update block properties
  block->UpdateProperties();
          
  // Insert block into database
  database->InsertBlock(block);
            
  // Insert block into node
  node->InsertBlock(block);
          
  // Update node properties
  node->UpdateProperties();
          
  // Release block
  database->ReleaseBlock(block);

  // Return node
  return node;
}
#endif



static R3SurfelNode *
LoadSurfels(R3SurfelScene *scene, const char *surfels_filename, 
  const char *object_name, const char *parent_object_name,
  const char *node_name, const char *parent_node_name)
{
  // Start statistics
  RNTime start_time;
  start_time.Read();

  // Get surfel tree
  R3SurfelTree *tree = scene->Tree();
  if (!tree) {
    RNFail("Scene has no tree\n");
    return NULL;
  }    

  // Get surfel database
  R3SurfelDatabase *database = tree->Database();
  if (!database) {
    RNFail("Scene has no database\n");
    return NULL;
  }    

  // Find parent object
  R3SurfelObject *parent_object = NULL;
  if (parent_object_name && 
      strcmp(parent_object_name, "None") && 
      strcmp(parent_object_name, "none") && 
      strcmp(parent_object_name, "NONE")) {
    parent_object = scene->FindObjectByName(parent_object_name);
    if (!parent_object) {
      RNFail("Unable to find parent object with name %s\n", parent_object_name);
      return NULL;
    }
  }

  // Find parent node
  R3SurfelNode *parent_node = tree->FindNodeByName(parent_node_name);
  if (!parent_node) {
    RNFail("Unable to find parent node with name %s\n", parent_node_name);
    return NULL;
  }

  // Create node
  R3SurfelNode *node = new R3SurfelNode(node_name);
  if (!node) {
    RNFail("Unable to allocate node for %s\n", surfels_filename);
    return NULL;
  }

  // Insert node into tree
  tree->InsertNode(node, parent_node);

  // Create block
  R3SurfelBlock *block = new R3SurfelBlock();
  if (!block) {
    RNFail("Unable to allocate block for %s\n", surfels_filename);
    return NULL;
  }

  // Read block
  if (!block->ReadFile(surfels_filename)) {
    RNFail("Unable to read block from %s\n", surfels_filename);
    return NULL;
  }

  // Extract subset of surfels
  if (aerial_only || terrestrial_only) {
    // Create subset
    R3SurfelPointSet *subset = new R3SurfelPointSet();

    // Fill subset
    for (int i = 0; i < block->NSurfels(); i++) {
      const R3Surfel *surfel = block->Surfel(i);
      if (surfel->IsAerial() && terrestrial_only) continue;
      if (!surfel->IsAerial() && aerial_only) continue;
      R3SurfelPoint point(block, surfel);
      subset->InsertPoint(point);
    }

    // Replace block with subset
    // Note: it is important to delete subset first, since it references block
    R3SurfelBlock *subset_block = new R3SurfelBlock(subset);
    delete subset;
    delete block;
    block = subset_block;
  }

  // Update block properties
  block->UpdateProperties();

  // Insert block into database
  database->InsertBlock(block);

  // Insert block into node
  node->InsertBlock(block);

  // Update node properties
  node->UpdateProperties();

  // Create object
  if (object_name && 
      strcmp(object_name, "None") && 
      strcmp(object_name, "none") && 
      strcmp(object_name, "NONE") && 
      parent_object) {
    // Create object
    R3SurfelObject *object = new R3SurfelObject(object_name);
    if (!object) {
      RNFail("Unable to create object\n");
      return NULL;
    }

    // Insert object into scene
    scene->InsertObject(object, parent_object);

    // Insert node into object
    object->InsertNode(node);
      
    // Update object properties
    object->UpdateProperties();
  }

  // Release block
  database->ReleaseBlock(block);

  // Print statistics
  if (print_verbose) {
    printf("Loaded surfels from %s ...\n", surfels_filename);
    printf("  Time = %.2f seconds\n", start_time.Elapsed());
    printf("  # Surfels = %g\n", node->Complexity());
    fflush(stdout);
  }

  // Return node
  return node;
}



static int
LoadSurfelsList(R3SurfelScene *scene, const char *list_filename, 
  const char *parent_object_name, const char *parent_node_name)
{
  // Start statistics
  RNTime start_time;
  start_time.Read();
  int saved_print_verbose = print_verbose;
  print_verbose = print_debug;

  // Open file
  FILE *fp = fopen(list_filename, "r");
  if (!fp) {
    RNFail("Unable to open %s\n", list_filename);
    return 0;
  }

  // Read objects/nodes from file with list
  int count = 0;
  char buffer[4096];
  char node_filename[512];
  while (fgets(buffer, 4096, fp)) {
    if (buffer[0] == '#') continue;
    if (sscanf(buffer, "%s", node_filename) == (unsigned int) 1) {
      char *start = strrchr(node_filename, '/');
      start = (start) ? start+1 : node_filename;
      char node_name[1024];
      strncpy(node_name, start, 1024);
      char *end = strrchr(node_name, '.');
      if (end) *end = '\0';
      if (!LoadSurfels(scene, node_filename, node_name, parent_object_name, node_name, parent_node_name)) return 0;
      count++;
    }
  }

  // Close file
  fclose(fp);

  // Print statistics
  print_verbose = saved_print_verbose;
  if (print_verbose) {
    printf("Loaded surfels from %s ...\n", list_filename);
    printf("  Time = %.2f seconds\n", start_time.Elapsed());
    printf("  # Files = %d\n", count);
    fflush(stdout);
  }

  // Return success
  return 1;
}



static int
LoadSurfelsFromMesh(R3SurfelScene *scene, const char *mesh_filename, 
  const char *parent_object_name, const char *parent_node_name,
  RNLength surfel_spacing = 0.01)
{
  // Start statistics
  RNTime start_time;
  start_time.Read();
  int surfel_count = 0;
  int node_count = 0;

  // Get surfel tree
  R3SurfelTree *tree = scene->Tree();
  if (!tree) {
    RNFail("Scene has no tree\n");
    return 0;
  }    

  // Get surfel database
  R3SurfelDatabase *database = tree->Database();
  if (!database) {
    RNFail("Scene has no database\n");
    return 0;
  }    

  // Find parent object
  R3SurfelObject *parent_object = scene->FindObjectByName(parent_object_name);
  if (!parent_object) {
    RNFail("Unable to find parent object with name %s\n", parent_object_name);
    return 0;
  }

  // Find parent node
  R3SurfelNode *parent_node = tree->FindNodeByName(parent_node_name);
  if (!parent_node) {
    RNFail("Unable to find parent node with name %s\n", parent_node_name);
    return 0;
  }

  // Read mesh file
  R3Mesh mesh;
  if (!mesh.ReadFile(mesh_filename)) {
    RNFail("Unable to read mesh from %s\n", mesh_filename);
    return 0;
  }

  // Create surfels for every mesh segment
  RNSeedRandomScalar();
  RNArray<R3SurfelLabel *> labels;
  RNArray<RNArray<R3Surfel *> *> surfels;
  R3Point mesh_centroid = mesh.Centroid();
  for (int i = 0; i < mesh.NFaces(); i++) {
    R3MeshFace *face = mesh.Face(i);
    int segment_index = mesh.FaceSegment(face);
    if (segment_index < 0) segment_index = mesh.FaceMaterial(face);
    if (segment_index < 0) segment_index = 0;

    // Get label
    int label_index = mesh.FaceCategory(face);
    if ((label_index >= 0) && (label_index < scene->NLabels())) {
      R3SurfelLabel *label = scene->Label(label_index);
      labels[segment_index] = label;
    }

    // Get/create array of surfels for segment
    while (labels.NEntries() <= segment_index) labels.Insert(NULL);
    while (surfels.NEntries() <= segment_index) surfels.Insert(NULL);
    RNArray<R3Surfel *> *segment_surfels = surfels[segment_index];
    if (!segment_surfels) {
      segment_surfels = new RNArray<R3Surfel *>();
      surfels[segment_index] = segment_surfels;
    }

    // Get vertex positions
    R3MeshVertex *v0 = mesh.VertexOnFace(face, 0);
    R3MeshVertex *v1 = mesh.VertexOnFace(face, 1);
    R3MeshVertex *v2 = mesh.VertexOnFace(face, 2);
    R3Point p0 = mesh.VertexPosition(v0);
    R3Point p1 = mesh.VertexPosition(v1);
    R3Point p2 = mesh.VertexPosition(v2);
    // const R3Vector& n0 = mesh.VertexNormal(v0);
    // const R3Vector& n1 = mesh.VertexNormal(v1);
    // const R3Vector& n2 = mesh.VertexNormal(v2);
    const RNRgb& c0 = mesh.VertexColor(v0);
    const RNRgb& c1 = mesh.VertexColor(v1);
    const RNRgb& c2 = mesh.VertexColor(v2);

    // Translate vertices by mesh centroid
    p0 -= mesh_centroid.Vector();
    p1 -= mesh_centroid.Vector();
    p2 -= mesh_centroid.Vector();

    // Determine number of surfels for face
    RNScalar surfels_per_area = 4.0 / (RN_PI * surfel_spacing * surfel_spacing);
    RNScalar ideal_face_nsurfels = surfels_per_area * mesh.FaceArea(face);
    int face_nsurfels = (int) ideal_face_nsurfels;
    RNScalar remainder = ideal_face_nsurfels - face_nsurfels;
    if (remainder > RNRandomScalar()) face_nsurfels++;

    // Generate random surfels on face
    for (int j = 0; j < face_nsurfels; j++) {
      RNScalar r1 = sqrt(RNRandomScalar());
      RNScalar r2 = RNRandomScalar();
      RNScalar t0 = (1.0 - r1);
      RNScalar t1 = r1 * (1.0 - r2);
      RNScalar t2 = r1 * r2;
      R3Point position = t0*p0 + t1*p1 + t2*p2;
      R3Vector normal = mesh.FaceNormal(face); // t0*n0 + t1*n1 + t2*n2; normal.Normalize();
      RNRgb color = t0*c0 + t1*c1 + t2*c2; 
      R3Surfel *surfel = new R3Surfel(position.X(), position.Y(), position.Z(),
        normal.X(), normal.Y(), normal.Z(), surfel_spacing,
        255*color.R(), 255*color.G(), 255*color.B(), R3_SURFEL_AERIAL_FLAG);
      segment_surfels->Insert(surfel);
    }
  }

  // Create surfel blocks/nodes
  for (int i = 0; i < surfels.NEntries(); i++) {
    R3SurfelLabel *label = labels.Kth(i);
    RNArray<R3Surfel *> *segment_surfels = surfels.Kth(i);
    if (!segment_surfels) continue;

    // Compute segment name
    char segment_name[1024];
    sprintf(segment_name, "MESH_SEGMENT_%d", i);

    // Create object
    R3SurfelObject *object = new R3SurfelObject(segment_name);
     if (!object) {
      RNFail("Unable to allocate object\n");
      return 0;
    }
   
    // Insert object into scene
    scene->InsertObject(object, parent_object);

    // Create node
    R3SurfelNode *node = new R3SurfelNode(segment_name);
    if (!node) {
      RNFail("Unable to allocate node\n");
      return 0;
    }
            
    // Insert node into tree
    tree->InsertNode(node, parent_node);
          
    // Insert node into object
    object->InsertNode(node);

    // Create block
    RNArray<const R3Surfel *>& tmp = *((RNArray<const R3Surfel *> *) segment_surfels);
    R3SurfelBlock *block = new R3SurfelBlock(tmp, mesh_centroid);
    if (!block) {
      RNFail("Unable to allocate block\n");
      return 0;
    }
    
    // Update block properties
    block->UpdateProperties();
          
    // Insert block into database
    database->InsertBlock(block);
            
    // Insert block into node
    node->InsertBlock(block);
          
    // Update node properties
    node->UpdateProperties();
          
    // Update object properties
    object->UpdateProperties();

    // Create label assignment
    if (label) {
      int originator = R3_SURFEL_LABEL_ASSIGNMENT_HUMAN_ORIGINATOR;
      R3SurfelLabelAssignment *assignment = new R3SurfelLabelAssignment(object, label, 1.0, originator);
      scene->InsertLabelAssignment(assignment);
    }

    // Release block
    database->ReleaseBlock(block);

    // Delete segment surfel arrays    
    for (int j = 0; j < segment_surfels->NEntries(); j++)
      delete segment_surfels->Kth(j);
    delete segment_surfels;

    // Increment statistics
    surfel_count += block->NSurfels();
    node_count++;
  }
  
  // Print statistics
  if (print_verbose) {
    printf("Loaded surfels from mesh %s ...\n", mesh_filename);
    printf("  Time = %.2f seconds\n", start_time.Elapsed());
    printf("  # Objects = %d\n", node_count);
    printf("  # Surfels = %d\n", surfel_count);
    fflush(stdout);
  }

  // Return success
  return 1;
}



#if 0
int MPHouse::
ReadCategoryFile(R3SurfelScene *scene, const char *filename, const char *root_name)
{
  // Start statistics
  RNTime start_time;
  start_time.Read();
  int count = 0;

  // Find root
  R3SurfelLabel *root = NULL;
  if (root_name && strcmp(root_name, "Null")) {
    root = scene->FindLabelByName(root_name);
    if (!root) {
      RNFail("Unable to find root label %s\n", root_name);
      return 0;
    }
  }
 
  // Open file
  FILE* fp = fopen(filename, "r");
  if (!fp) {
    RNFail("Unable to open category file %s\n", filename);
    return 0;
  }

  // Read keys from first line
  RNArray<char *> keys;
  char key_buffer[4096];
  if (fgets(key_buffer, 4096, fp)) {
    char *token = key_buffer;
    while (*token &&  (*token != '\r') && (*token != '\n')) {
      while (*token == ' ') token++;
      keys.Insert(token);
      while (*token && (*token != '\t') && (*token != '\r') && (*token != '\n')) token++;
      *(token++) = '\0';
    }
  }
  
  // Extract indices of interesting info
  int label_id_k = -1;
  int label_name_k = -1;
  int mpcat40_id_k = -1;
  int mpcat40_name_k = -1;
  for (int i = 0; i < keys.NEntries(); i++) {
    if (!strcmp(keys[i], "index")) label_id_k = i;
    else if (!strcmp(keys[i], "raw_category")) label_name_k = i; 
    else if (!strcmp(keys[i], "mpcat40index")) mpcat40_id_k = i; 
    else if (!strcmp(keys[i], "mpcat40")) mpcat40_name_k = i; 
  }

  // Check if found id field in header
  if ((label_id_k < 0) && (mpcat40_id_k < 0) {
    RNFail("Did not find index or mpcat40index in header of %s\n", filename);
    return 0;
  }

  // Check if found name field in header
  if ((label_name_k < 0) && (mpcat40_name_k < 0) {
    RNFail("Did not find raw_category or mpcat40 in header of %s\n", filename);
    return 0;
  }

  // Read subsequent lines of file
  char value_buffer[4096];
  while (fgets(value_buffer, 4096, fp)) {
    // Read values
    RNArray<char *> values;
    char *token = value_buffer;
    while (*token &&  (*token != '\r') && (*token != '\n')) {
      while (*token == ' ') token++;
      values.Insert(token);
      while (*token && (*token != '\t') && (*token != '\r') && (*token != '\n')) token++;
      *(token++) = '\0';
    }

    // Create category
    int label_id = (mpcat40_name_k >= 0) ? atoi(values[label_id_k]) : -1;
    char *label_name =(mpcat40_name_k >= 0) ? RNStrdup(values[label_name_k]) : NULL;
    int mpcat40_id = (mpcat40_name_k >= 0) ? atoi(values[mpcat40_id_k]) : -1;
    char *mpcat40_name = (mpcat40_name_k >= 0) ? RNStrdup(values[mpcat40_name_k]) : NULL;
    int id = (mpcat40_id > 0) ? mpcat40_id : label_id;
    char *name = (mpcat40_name) ? mpcat40_name : label_name;
              
    // Insert labels up to id, so that id matches index
    while (scene->NLabels() < id) {
      scene->InsertLabel(new R3SurfelLabel());
    }
    
    // Set label info
    R3SurfelLabel *label = scene->Label(id);
    label->SetName(name);
    label->SetIdentifier(id);
    label->SetColor(RNRgb(RNRandomScalar(), RNRandomScalar(), RNRandomScalar()));
    
    // Update stats
    count++;
  }
  
  // Close file
  fclose(fp);

  // Print statistics
  if (print_verbose) {
    printf("Read categories from %s ...\n", filename);
    printf("  Time = %.2f seconds\n", start_time.Elapsed());
    printf("  # Labels = %d\n", count);
    fflush(stdout);
  }

  // Return success
  return 1;
}
#endif



static int
LoadLabelList(R3SurfelScene *scene, const char *list_filename, const char *root_name)
{
  // Start statistics
  RNTime start_time;
  start_time.Read();

  // Open file
  FILE *fp = fopen(list_filename, "r");
  if (!fp) {
    RNFail("Unable to open %s\n", list_filename);
    return 0;
  }

  // Find root
  R3SurfelLabel *root = NULL;
  if (root_name && strcmp(root_name, "Null")) {
    root = scene->FindLabelByName(root_name);
    if (!root) {
      RNFail("Unable to find root label %s\n", root_name);
      return 0;
    }
  }
 
  // Read labels from file with list
  int count = 0;
  double r, g, b;
  int identifier, visibility;
  char assignment_keystroke[64];
  char label_name[4096], parent_name[4096], buffer[16384];
  while (fgets(buffer, 16384, fp)) {
    char *bufferp = buffer;
    while (*bufferp && isspace(*bufferp)) bufferp++;
    if (*bufferp == '\0') continue;
    if (*bufferp == '#') continue;
    if (sscanf(buffer, "%s%d%s%s%d%lf%lf%lf", label_name, &identifier, assignment_keystroke, parent_name, &visibility, &r, &g, &b) != (unsigned int) 8) {
      RNFail("Invalid format for label %d in %s\n", count, list_filename);
      return 0;
    }
          
    // Check if label already exists
    if (scene->FindLabelByName(label_name)) continue;

    // Create label
    R3SurfelLabel *label = new R3SurfelLabel(label_name);
    if (assignment_keystroke[0] != '-') label->SetAssignmentKeystroke(assignment_keystroke[0]);
    label->SetIdentifier(identifier);
    label->SetColor(RNRgb(r, g, b));
    
    // Find parent
    R3SurfelLabel *parent = NULL;
    if (!strcmp(parent_name, "Null")) parent = root;
    else {
      parent = scene->FindLabelByName(parent_name);
      if (!parent) {
        RNFail("Unable to find label's parent (%s) in label %d of %s\n", parent_name, count, list_filename);
        return 0;
      }
    }
    
    // Insert into scene
    scene->InsertLabel(label, parent);
    
    // Update stats
    count++;
  }

  // Close file
  fclose(fp);

  // Print statistics
  if (print_verbose) {
    printf("Loaded labels from %s ...\n", list_filename);
    printf("  Time = %.2f seconds\n", start_time.Elapsed());
    printf("  # Labels = %d\n", count);
    fflush(stdout);
  }

  // Return success
  return 1;
}



static int
LoadAssignmentList(R3SurfelScene *scene, const char *list_filename)
{
  // Start statistics
  RNTime start_time;
  start_time.Read();
  int count = 0;

  // Open file
  FILE *fp = fopen(list_filename, "r");
  if (!fp) {
    RNFail("Unable to open %s\n", list_filename);
    return 0;
  }

  // Read labels from file with list
  char buffer[4096];
  char object_name[4096];
  char label_name[4096];
  char originator_str[4096];
  double confidence;
  while (fgets(buffer, 4096, fp)) {
    if (buffer[0] == '\0') continue;
    if (buffer[0] == '#') continue;
    if (sscanf(buffer, "%s%s%lf%s", object_name, label_name, &confidence, originator_str) == 4) {
      // Find object
      R3SurfelObject *object = scene->FindObjectByName(object_name);
      if (!object) {
        RNFail("Unable to find object %s in assignments file %s\n", object_name, list_filename);
        return 0;
      }

      // Find label
      R3SurfelLabel *label = scene->FindLabelByName(label_name);
      if (!label) {
        RNFail("Unable to find label %s in assignments file %s\n", label_name, list_filename);
        return 0;
      }

      // Create assignment
      int originator = R3_SURFEL_LABEL_ASSIGNMENT_MACHINE_ORIGINATOR;
      if (!strcmp(originator_str, "Human")) originator = R3_SURFEL_LABEL_ASSIGNMENT_HUMAN_ORIGINATOR;
      else if (!strcmp(originator_str, "GroundTruth")) originator = R3_SURFEL_LABEL_ASSIGNMENT_GROUND_TRUTH_ORIGINATOR;
      R3SurfelLabelAssignment *assignment = new R3SurfelLabelAssignment(object, label, confidence, originator);
      scene->InsertLabelAssignment(assignment);
      count++;
    }
  }

  // Close file
  fclose(fp);

  // Print statistics
  if (print_verbose) {
    printf("Loaded assignments from %s ...\n", list_filename);
    printf("  Time = %.2f seconds\n", start_time.Elapsed());
    printf("  # Assignments = %d\n", count);
    fflush(stdout);
  }

  // Return success
  return 1;
}



static int
LoadFeatureList(R3SurfelScene *scene, const char *list_filename)
{
  // Start statistics
  RNTime start_time;
  start_time.Read();
  int count = 0;

  // Open file
  FILE *fp = fopen(list_filename, "r");
  if (!fp) {
    RNFail("Unable to open %s\n", list_filename);
    return 0;
  }

  // Read features from file with list
  char buffer[4096];
  char type[1024], name[1024], filename[1024];
  double minimum, maximum, weight;
  while (fgets(buffer, 4096, fp)) {
    if (buffer[0] == '\0') continue;
    if (buffer[0] == '#') continue;
    if (sscanf(buffer, "%s%s%lf%lf%lf%s", type, name, &minimum, &maximum, &weight, filename) == 6) {
      // Create feature of appropriate type
      if (!strcmp(type, "PointSet")) {
        R3SurfelPointSetFeature *feature = new R3SurfelPointSetFeature(name, minimum, maximum, weight);
        scene->InsertFeature(feature);
      }
      else if (!strcmp(type, "OverheadGrid")) {
        R3SurfelOverheadGridFeature *feature = new R3SurfelOverheadGridFeature(filename, name, minimum, maximum, weight);
        scene->InsertFeature(feature);
      }
      else {
        R3SurfelFeature *feature = new R3SurfelFeature(name, minimum, maximum, weight);
        scene->InsertFeature(feature);
      }
      count++;
    }
  }

  // Close file
  fclose(fp);

  // Print statistics
  if (print_verbose) {
    printf("Loaded features from %s ...\n", list_filename);
    printf("  Time = %.2f seconds\n", start_time.Elapsed());
    printf("  # Features = %d\n", count);
    fflush(stdout);
  }

  // Return success
  return 1;
}



static int
LoadScene(R3SurfelScene *scene1, 
  const char *scene_filename, const char *database_filename, 
  const char *parent_object_name, const char *parent_label_name, const char *parent_node_name)
{
  // Start statistics
  RNTime start_time;
  start_time.Read();

  // Find parent surfel tree node in scene1
  R3SurfelNode *parent_node = scene1->Tree()->FindNodeByName(parent_node_name);
  if (!parent_node) {
    RNFail("Unable to find parent node with name %s\n", parent_node_name);
    return 0;
  }

  // Find parent object in scene1
  R3SurfelObject *parent_object = scene1->FindObjectByName(parent_object_name);
  if (!parent_object) {
    RNFail("Unable to find parent object with name %s\n", parent_object_name);
    return 0;
  }

  // Find parent label in scene1
  R3SurfelLabel *parent_label = scene1->FindLabelByName(parent_label_name);
  if (!parent_label) {
    RNFail("Unable to find parent label with name %s\n", parent_label_name);
    return 0;
  }

  // Allocate scene2
  R3SurfelScene *scene2 = new R3SurfelScene();
  if (!scene2) {
    RNFail("Unable to allocate scene\n");
    return 0;
  }

  // Open scene2
  if (!scene2->OpenFile(scene_filename, database_filename, "r", "r")) {
    delete scene2;
    return 0;
  }

  // Insert scene2 into scene1
  scene1->InsertScene(*scene2, parent_object, parent_label, parent_node);

  // Print statistics
  if (print_verbose) {
    printf("Loaded scene from %s and %s ...\n", scene_filename, database_filename);
    printf("  Time = %.2f seconds\n", start_time.Elapsed());
    printf("  # Objects = %d\n", scene2->NObjects());
    printf("  # Labels = %d\n", scene2->NLabels());
    printf("  # Assignments = %d\n", scene2->NLabelAssignments());
    printf("  # Nodes = %d\n", scene2->Tree()->NNodes());
    printf("  # Blocks = %d\n", scene2->Tree()->Database()->NBlocks());
    printf("  # Surfels = %lld\n", scene2->Tree()->Database()->NSurfels());
    fflush(stdout);
  }

  // Close scene
  if (!scene2->CloseFile()) return 0;

  // Delete scene
  delete scene2;

  // Return success
  return 1;
}



////////////////////////////////////////////////////////////////////////
// ALIGNMENT FUNCTIONS
////////////////////////////////////////////////////////////////////////

static int
TransformWithConfigurationFile(R3SurfelScene *scene, const char *filename, RNBoolean invert = FALSE)
{
  // Start statistics
  RNTime start_time;
  start_time.Read();
  int count = 0;
  
  // Get surfel tree
  R3SurfelTree *tree = scene->Tree();
  if (!tree) {
    RNFail("Scene has no surfel tree\n");
    return 0;
  }    

  // Open file
  FILE *fp = fopen(filename, "r");
  if (!fp) { 
    RNFail("Unable to open extrinsics file %s\n", filename); 
    return 0; 
  }

  // Parse file
  char buffer[4096];
  int line_number = 0;
  while (fgets(buffer, 4096, fp)) {
    char cmd[4096];
    line_number++;
    if (sscanf(buffer, "%s", cmd) != (unsigned int) 1) continue;
    if (cmd[0] == '#') continue;

    // Check cmd
    if (!strcmp(cmd, "scan") || !strcmp(cmd, "image") || !strcmp(cmd, "frame")) {
      // Parse image name and alignment transformation
      RNScalar m[16], depth_timestamp, color_timestamp;
      char depth_name[1024], color_name[1024];
      if (!strcmp(cmd, "frame")) {
        if (sscanf(buffer, "%s%s%lf%s%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf", cmd, 
           depth_name, &depth_timestamp, color_name, &color_timestamp,
           &m[0], &m[1], &m[2], &m[3], &m[4], &m[5], &m[6], &m[7], 
           &m[8], &m[9], &m[10], &m[11], &m[12], &m[13], &m[14], &m[15]) != (unsigned int) 21) {
          RNFail("Error parsing line %d of %s\n", line_number, filename);
          return 0;
        }
      }
      else {
        if (sscanf(buffer, "%s%s%s%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf", cmd, 
           depth_name, color_name,
           &m[0], &m[1], &m[2], &m[3], &m[4], &m[5], &m[6], &m[7], 
           &m[8], &m[9], &m[10], &m[11], &m[12], &m[13], &m[14], &m[15]) != (unsigned int) 19) {
          RNFail("Error parsing line %d of %s\n", line_number, filename);
          return 0;
        }
      }

      // Get transformation
      R3Affine transformation(R4Matrix(m), 0);
      if (invert) transformation.Invert();

      // Get image name
      char tmp[2048];
      strncpy(tmp, depth_name, 1048);
      char *image_name = strrchr(tmp, '/');
      if (!image_name) image_name = tmp;
      char *endp = strrchr(image_name, '.');
      if (endp) *endp = '\0';
      char node_name[4096];
      sprintf(node_name, "SCAN:%s", image_name);

      // Find node
      R3SurfelNode *node = tree->FindNodeByName(node_name);
      R3SurfelScan *scan = scene->FindScanByName(image_name);
      R3SurfelImage *image = scene->FindImageByName(image_name);
      
      // Transform image
      if (image) {
        R3CoordSystem pose = image->Pose();
        pose.Transform(transformation);
        image->SetPose(pose);
      }
      
      // Transform scan
      if (scan) {
        R3CoordSystem pose = scan->Pose();
        pose.Transform(transformation);
        scan->SetPose(pose);
      }
      
      // Transform node and all its decendents
      if (node) {
        RNArray<R3SurfelNode *> stack;
        stack.Insert(node);
        while (!stack.IsEmpty()) {
          R3SurfelNode *node = stack.Tail();
          stack.RemoveTail();
          node->Transform(transformation);
          for (int i = 0; i < node->NParts(); i++) {
            stack.InsertTail(node->Part(i));
          }
        }
      }

      // Update statistics
      count++;
    }
  }

  // Close configuration file
  fclose(fp);

  // Print statistics
  if (print_verbose) {
    printf("Tranformed nodes ...\n");
    printf("  Time = %.2f seconds\n", start_time.Elapsed());
    printf("  # Nodes = %d\n", count);
    fflush(stdout);
  }

  // Return success
  return 1;
}



////////////////////////////////////////////////////////////////////////
// TRANSFER FUNCTIONS
////////////////////////////////////////////////////////////////////////

#if 0

static int
TransferLabels(R3SurfelScene *scene1, R3SurfelScene *scene2)
{
  // Copy labels from scene2 to scene1
  for (int i = 0; i < scene2->NLabels(); i++) {
    R3SurfelLabel *label2 = scene2->Label(i);
    R3SurfelLabel *label1 = scene1->FindLabelByName(label2->Name());
    if (!label1) scene1->InsertLabel(new R3SurfelLabel(*label2));
  }

  // Return success
  return 1;
}



static int
TransferObjects(R3SurfelScene *scene1, R3SurfelScene *scene2)
{
  // Initialize an object index

  // Copy objects from scene2 to scene1
  R3SurfelObject *root = scene1->RootObject();
  for (int i = 0; i < scene2->NObjects(); i++) {
    R3SurfelObject *object2 = scene2->Object(i);
    if (object2 == scene2->RootObject()) continue;
    R2SurfelObject *object1 = new R3SurfelObject(object2->Name());
    scene1->InsertObject(object1, root);
    // Extract surfels nearby ones in object2 and move them into object1
  }

  // Copy object hierarchy from scene2 to scene1
  // Must be done in separate pass in case parents are created after children
  for (int i = 0; i < scene2->NObjects(); i++) {
    R3SurfelObject *object2 = scene2->Object(i);
    R3SurfelObject *object1 = objects_index[i];
    R3SurfelObject *parent2 = object2->Parent();
    R3SurfelObject *parent1 = objects_index[parent2->SceneIndex()];
    object1->SetParent(parent1);
  }

  // Copy label assignments from scene2 to scene1
  for (int i = 0; i < scene2->NLabelAssignments(); i++) {
    R3SurfelLabelAssignment *assignment2 = scene2->LabelAssignment(i);
    R3SurfelObject *object2 = assignment2->Object();
    R3SurfelObject *object1 = objects_index[object2->SceneIndex()];
    R3SurfelLabel *label2 = assignment2->Label();
    R3SurfelLabel *label1 = scene1->FindLabelByName(label2->Name());
    R3SurfelLabelAssignment *assignment1 = new R3SurfelLabelAssignment(object1, label1, assignment2->Confidence(), assignment2->Originator());
    scene1->InsertLabelAssignment(assignment1);
  }

  // Return success
  return 1;
}



static int
TransferLabels(R3SurfelScene *scene1,  
  const char *label_scene_filename, const char *label_database_filename)
{
  // Start statistics
  RNTime start_time;
  start_time.Read();

  // Allocate scene2
  R3SurfelScene *scene2 = new R3SurfelScene();
  if (!scene2) {
    RNFail("Unable to allocate scene\n");
    return 0;
  }

  // Open scene2
  if (!scene2->OpenFile(label_scene_filename, label_database_filename, "r", "r")) {
    delete scene2;
    return 0;
  }

  // Transfer labels
  if (!TransferLabels(scene1, scene2)) return 0;

  // Transfer objects
  if (!TransferObjects(scene1, scene2)) return 0;

  // Print statistics
  if (print_verbose) {
    printf("Transferred labels from %s and %s ...\n", scene_filename, database_filename);
    printf("  Time = %.2f seconds\n", start_time.Elapsed());
    printf("  # Objects = %d\n", scene2->NObjects());
    printf("  # Labels = %d\n", scene2->NLabels());
    printf("  # Assignments = %d\n", scene2->NLabelAssignments());
    printf("  # Nodes = %d\n", scene2->Tree()->NNodes());
    printf("  # Blocks = %d\n", scene2->Tree()->Database()->NBlocks());
    printf("  # Surfels = %lld\n", scene2->Tree()->Database()->NSurfels());
    fflush(stdout);
  }

  // Close scene
  if (!scene2->CloseFile()) return 0;

  // Delete scene
  delete scene2;

  // Return success
  return 1;
}

#endif



////////////////////////////////////////////////////////////////////////
// MASKING FUNCTIONS
////////////////////////////////////////////////////////////////////////

static int
Mask(R3SurfelScene *scene, const char *node_name, R3SurfelConstraint *constraint)
{
  // Start statistics
  RNTime start_time;
  start_time.Read();

  // Get useful variables
  R3SurfelTree *tree = scene->Tree();
  if (!tree) {
    RNFail("Scene has no surfel tree\n");
    return 0;
  }    

  // Get surfel database
  R3SurfelDatabase *database = tree->Database();
  if (!database) {
    RNFail("Tree has no surfel database\n");
    return 0;
  }    

  // Find node
  R3SurfelNode *node = tree->RootNode();
  if (strcmp(node_name, "All") && strcmp(node_name, "Root")) {
    node = tree->FindNodeByName(node_name);
    if (!node) {
      RNFail("Unable to find node with name %s\n", node_name);
      return 0;
    }
  }

  // Remove surfels not satistfying constraint
  RNArray<R3SurfelNode *> remove_nodes;
  tree->SplitLeafNodes(node, *constraint, NULL, &remove_nodes);
  for (int i = 0; i < remove_nodes.NEntries(); i++) {
    R3SurfelNode *node = remove_nodes.Kth(i);

    // Make array of blocks
    RNArray<R3SurfelBlock *> blocks;
    while (node->NBlocks() > 0) {
      R3SurfelBlock *block = node->Block(0);
      node->RemoveBlock(block);
      blocks.Insert(block);
    }

    // Remove/delete node
    tree->RemoveNode(node);
    delete node;

    // Remove/delete blocks
    for (int j = 0; j < blocks.NEntries(); j++) {
      R3SurfelBlock *block = blocks.Kth(j);
      database->RemoveBlock(block);
      delete block;
    }
  }

  // Return success
  return 1;
}



////////////////////////////////////////////////////////////////////////
// MULTIRESOLUTION FUNCTIONS
////////////////////////////////////////////////////////////////////////

static int
SplitSurfelTreeNodes(R3SurfelScene *scene, const char *node_name,
  int max_parts_per_node, int max_blocks_per_node, 
  RNScalar max_node_complexity, RNScalar max_block_complexity,
  RNLength max_leaf_extent, RNLength max_block_extent, 
  int max_levels)
{
  // Start statistics
  RNTime start_time;
  start_time.Read();

  // Get surfel tree
  R3SurfelTree *tree = scene->Tree();
  if (!tree) {
    RNFail("Scene has no surfel tree\n");
    return 0;
  }    

  // Find node
  R3SurfelNode *node = NULL;
  if (strcmp(node_name, "All")) {
    node = tree->FindNodeByName(node_name);
    if (!node) {
      RNFail("Unable to find node with name %s\n", node_name);
      return 0;
    }
  }

  // Check node
  if (node) {
    // Split nodes
    tree->SplitNodes(node,
      max_parts_per_node, max_blocks_per_node, 
      max_node_complexity, max_block_complexity, 
      max_leaf_extent, max_block_extent, 
      max_levels);

    // Print statistics
    if (print_verbose) {
      printf("Split nodes starting at %s ...\n", node_name);
      printf("  Time = %.2f seconds\n", start_time.Elapsed());
      printf("  # Nodes = %d\n", tree->NNodes());
      printf("  # Blocks = %d\n", tree->Database()->NBlocks());
      fflush(stdout);
    }
  }
  else {
    // Split nodes
    tree->SplitNodes(
      max_parts_per_node, max_blocks_per_node, 
      max_node_complexity, max_block_complexity, 
      max_leaf_extent, max_block_extent, 
      max_levels);

    // Print statistics
    if (print_verbose) {
      printf("Split all nodes  ...\n");
      printf("  Time = %.2f seconds\n", start_time.Elapsed());
      printf("  # Nodes = %d\n", tree->NNodes());
      printf("  # Blocks = %d\n", tree->Database()->NBlocks());
      fflush(stdout);
    }
  }

  // Return success
  return 1;
}



static int
CreateMultiresolutionNodes(R3SurfelScene *scene, const char *node_name,
  RNScalar min_complexity, RNScalar min_resolution, RNScalar min_multiresolution_factor)
{
  // Start statistics
  RNTime start_time;
  start_time.Read();

  // Get surfel tree
  R3SurfelTree *tree = scene->Tree();
  if (!tree) {
    RNFail("Scene has no surfel tree\n");
    return 0;
  }    

  // Temporary
  if (strcmp(node_name, "All")) {
    RNFail("-create_multiresolution_nodes only supported for All nodes\n");
    return 0;
  }

  // Create multiresolution nodes
  tree->CreateMultiresolutionNodes(min_complexity, min_resolution, min_multiresolution_factor);

  // Print statistics
  if (print_verbose) {
    printf("Created multiresolution nodes  ...\n");
    printf("  Time = %.2f seconds\n", start_time.Elapsed());
    printf("  # Nodes = %d\n", tree->NNodes());
    fflush(stdout);
  }

  // Return success
  return 1;
}



static int
CreateMultiresolutionBlocks(R3SurfelScene *scene, const char *node_name,
  RNScalar multiresolution_factor, RNScalar max_node_complexity)
{
  // Start statistics
  RNTime start_time;
  start_time.Read();

  // Get surfel tree
  R3SurfelTree *tree = scene->Tree();
  if (!tree) {
    RNFail("Scene has no surfel tree\n");
    return 0;
  }    

  // Find node
  R3SurfelNode *node = NULL;
  if (strcmp(node_name, "All")) {
    node = tree->FindNodeByName(node_name);
    if (!node) {
      RNFail("Unable to find node with name %s\n", node_name);
      return 0;
    }
  }

  // Check node
  if (node) {
    // Create multiresolution starting at node
    tree->CreateMultiresolutionBlocks(node, multiresolution_factor, max_node_complexity);

    // Print statistics
    if (print_verbose) {
      printf("Created multiresolution blocks for nodes starting at %s ...\n", node_name);
      printf("  Time = %.2f seconds\n", start_time.Elapsed());
      printf("  # Nodes = %d\n", tree->NNodes());
      printf("  # Blocks = %d\n", tree->Database()->NBlocks());
      fflush(stdout);
    }
  }
  else {
    // Create multiresolution starting at all root nodes
    tree->CreateMultiresolutionBlocks(multiresolution_factor, max_node_complexity);

    // Print statistics
    if (print_verbose) {
      printf("Created multiresolution blocks for all nodes  ...\n");
      printf("  Time = %.2f seconds\n", start_time.Elapsed());
      printf("  # Nodes = %d\n", tree->NNodes());
      printf("  # Blocks = %d\n", tree->Database()->NBlocks());
      fflush(stdout);
    }
  }

  // Return success
  return 1;
}



////////////////////////////////////////////////////////////////////////
// RELATIONSHIP FUNCTIONS
////////////////////////////////////////////////////////////////////////

#define USE_GRID



#ifndef USE_GRID

struct SurfelPointCompatibilityParameters {
  RNScalar max_gap_distance;
  RNScalar max_normal_angle;
  RNScalar max_plane_offset;
};


static int
AreSurfelPointsCompatible(R3SurfelPoint *point0, R3SurfelPoint *point1, void *data)
{
  // Get parameters
  SurfelPointCompatibilityParameters *p = (SurfelPointCompatibilityParameters *) data;
  if (!p) return 1;

  // Check normal angle
  if (p->max_normal_angle > 0) {
    const R3Vector& normal0 = point0->Normal();
    const R3Vector& normal1 = point1->Normal();
    RNAngle normal_angle = R3InteriorAngle(normal0, normal1);
    if (normal_angle > p->max_normal_angle) return 0;
  }

  // Passed all tests
  return 1;
}

#endif



static int
CreateOverlapObjectRelationships(R3SurfelScene *scene,
  RNLength max_gap_distance, RNLength max_plane_offset, RNAngle max_normal_angle,
  RNScalar min_overlap)
{
  // Check max gap
  if (max_gap_distance == 0) return 1;
  
  // Start statistics
  RNTime start_time;
  start_time.Read();
  int count = 0;
  if (print_verbose) {
    printf("Creating object overlap relationships ...\n");
    fflush(stdout);
  }

  // Create a pointset for each object
  R3SurfelPointSet **pointsets = new R3SurfelPointSet * [ scene->NObjects() ];
  RNScalar max_resolution = 8.0/(max_gap_distance*max_gap_distance);
  for (int i = 0; i < scene->NObjects(); i++) {
    R3SurfelObject *object = scene->Object(i);
    R3SurfelNodeSet nodes;
    pointsets[i] = new R3SurfelPointSet();
    nodes.InsertNodes(object, max_resolution);
    for (int j = 0; j < nodes.NNodes(); j++) {
      R3SurfelNode *node = nodes.Node(j);
      for (int k = 0; k < node->NBlocks(); k++) {
        R3SurfelBlock *block = node->Block(k);
        pointsets[i]->InsertPoints(block);
      }
    }
  }

  // Consider every object
  for (int i0 = 0; i0 < scene->NObjects(); i0++) {
    R3SurfelObject *object0 = scene->Object(i0);
    R3SurfelPointSet *pointset0 = pointsets[i0];
    if (pointset0->NPoints() == 0) continue;
    const R3Box& bbox0 = pointset0->BBox();

#ifdef USE_GRID
    // Create grid
    R3Grid *grid0 = CreateGrid(pointset0, max_gap_distance);
    if (!grid0) continue;
#else
    // Create kdtree
    RNArray<R3SurfelPoint *> array0;
    for (int j = 0; j < pointset0->NPoints(); j++) array0.Insert(pointset0->Point(j));
    R3Kdtree<R3SurfelPoint *> kdtree(array0, SurfelPointPosition);
#endif
                                     
    // Consider every other object
    for (int i1 = 0; i1 < scene->NObjects(); i1++) {
      R3SurfelObject *object1 = scene->Object(i1);
      if (object0 == object1) continue;
      R3SurfelPointSet *pointset1 = pointsets[i1];
      if (pointset1->NPoints() == 0) continue;
      const R3Box& bbox1 = pointset1->BBox();

      // Check bboxes for max_gap_distance
      if (R3Distance(bbox0, bbox1) > max_gap_distance) continue;

      // Check pointset1 for overlap with grid0
      int npoints = 0;
      for (int j1 = 0; j1 < pointset1->NPoints(); j1++) {
        R3SurfelPoint *point1 = pointset1->Point(j1);
#ifdef USE_GRID
        // Check if point1 is in marked grid cell
        const R3Point& world_position = point1->Position();
        R3Point grid_position = grid0->GridPosition(world_position);
        int ix = (int) (grid_position.X() + 0.5);
        if ((ix < 0) || (ix >= grid0->XResolution())) continue;
        int iy = (int) (grid_position.Y() + 0.5);
        if ((iy < 0) || (iy >= grid0->YResolution())) continue;
        int iz = (int) (grid_position.Z() + 0.5);
        if ((iz < 0) || (iz >= grid0->ZResolution())) continue;
        if (grid0->GridValue(ix, iy, iz) <= 0) continue;
        // grid0->SetGridValue(ix, iy, iz, 0);
#else
        // Check if point1 is near compatible point in kdtree
        SurfelPointCompatibilityParameters compatibility_parameters;
        compatibility_parameters.max_gap_distance = max_gap_distance;
        compatibility_parameters.max_normal_angle = max_normal_angle;
        compatibility_parameters.max_plane_offset = max_plane_offset;
        if (!kdtree.FindAny(point1, 0, max_gap_distance, AreSurfelPointsCompatible, &compatibility_parameters)) continue;
#endif
        npoints++;
      }

      // Compute/check overlap
      RNScalar overlap = (RNScalar) npoints / (RNScalar) pointset1->NPoints();

      // Create object relationship
      if (overlap > min_overlap) {
        RNScalar operands[] = { overlap };
        int noperands = sizeof(operands) / sizeof(RNScalar);
        R3SurfelObjectRelationship *relationship = new R3SurfelObjectRelationship(R3_SURFEL_OBJECT_OVERLAP_RELATIONSHIP, object0, object1, operands, noperands);
        scene->InsertObjectRelationship(relationship);
        count++;
      }
    }

#ifdef USE_GRID
    // Delete grid
    delete grid0;
#endif
  }

  // Delete pointsets
  for (int i = 0; i < scene->NObjects(); i++) {
    delete pointsets[i];
  }

  // Print statistics
  if (print_verbose) {
    printf("  Time = %.2f seconds\n", start_time.Elapsed());
    printf("  # Objects = %d\n", scene->NObjects());
    printf("  # Relationships = %d\n", count);
    fflush(stdout);
  }
  
  // Return success
  return 1;
}



#if 0
static int
CreatePlaneObjectRelationships(R3SurfelScene *scene, RNLength max_gap_distance, RNLength max_plane_offset, RNAngle max_normal_angle)
{
  // Consider all objects0
  for (int i0 = 0; i0 < scene->NObjects(); i0++) {
    R3SurfelObject *object0 = scene->Object(i0);
    R3SurfelObjectProperty *property0 = object0->FindObjectProperty(R3_SURFEL_OBJECT_PCA_PROPERTY);
    if (!property0) continue;
    R3Point centroid0(property->Operand(0), property->Operand(1), property->Operand(2));
    R3Vector majoraxis0(property->Operand(3), property->Operand(4), property->Operand(5));
    R3Vector minoraxis0(property->Operand(6), property->Operand(7), property->Operand(8));
    R3Vector normal0(property->Operand(9), property->Operand(10), property->Operand(11));
    R3Plane plane0(centroi0, normal0);
    R3Vector stddev0(property->Operand(12), property->Operand(13), property->Operand(14));
    R3Box extent0(property->Operand(15), property->Operand(16), property->Operand(17),
                  property->Operand(18), property->Operand(19), property->Operand(20));
    
    // Consider all other objects1
    for (int i1 = i0+1; i1 < scene->NObjects(); i1++) {
      R3SurfelObject *object1 = scene->Object(i1);
      R3SurfelObjectProperty *property1 = object1->FindObjectProperty(R3_SURFEL_OBJECT_PCA_PROPERTY);
      if (!property1) continue;
      R3Point centroid1(property->Operand(0), property->Operand(1), property->Operand(2));
      R3Vector majoraxis1(property->Operand(3), property->Operand(4), property->Operand(5));
      R3Vector minoraxis1(property->Operand(6), property->Operand(7), property->Operand(8));
      R3Vector normal1(property->Operand(9), property->Operand(10), property->Operand(11));
      R3Plane plane1(centroi1, normal1);
      R3Vector stddev1(property->Operand(12), property->Operand(13), property->Operand(14));
      R3Box extent1(property->Operand(15), property->Operand(16), property->Operand(17),
                    property->Operand(18), property->Operand(19), property->Operand(20));

      // Check if closest points are close enough 
      if (max_gap_distance > 0) {
        RNLength gap = xxx;
        if (gap > max_gap_distance) continue;
      }

      // Check if centroids are withing max offset of others' planes
      if (max_offset > 0) {
        RNLength offset0 = R3Distance(plane0, centroid1);
        if (offset0 > max_offset) continue;
        RNLength offset1 = R3Distance(plane1, centroid0);
        if (offset1 > max_offset) continue;
      }

      // Check if normals are within max angle
      if (max_angle) {
        RNAngle angle = R3InteriorAngle(normal0, normal1);
        if (angle > max_angle) continue;
      }
    }
  }
  
  // Return success
  return 1;
}
#endif



static int
CreateObjectRelationships(R3SurfelScene *scene,
  RNLength max_gap_distance, RNLength max_plane_offset, RNAngle max_normal_angle,
  RNScalar min_overlap)
{
  // Create object relationships
  if (!CreateOverlapObjectRelationships(scene, max_gap_distance, max_plane_offset, max_normal_angle, min_overlap)) return 0;

  // Return success
  return 1;
}



////////////////////////////////////////////////////////////////////////
// SEGMENTATION FUNCTIONS
////////////////////////////////////////////////////////////////////////

static int
CreateClusterObjects(R3SurfelScene *scene, 
  const char *parent_object_name, const char *parent_node_name, const char *source_node_name, 
  int max_neighbors, RNLength max_neighbor_distance, 
  RNLength max_offplane_distance, RNAngle max_normal_angle,
  int min_points_per_object, RNLength chunk_size)
{
  // Start statistics
  RNTime start_time;
  start_time.Read();
  if (print_verbose) {
    printf("Creating cluster objects ...\n");
    fflush(stdout);
  }

  // Find parent object
  R3SurfelObject *parent_object = scene->FindObjectByName(parent_object_name);
  if (!parent_object) {
    RNFail("Unable to find object with name %s\n", parent_object_name);
    return 0;
  }

  // Find parent node
  R3SurfelNode *parent_node = scene->Tree()->FindNodeByName(parent_node_name);
  if (!parent_node) {
    RNFail("Unable to find node with name %s\n", parent_node_name);
    return 0;
  }

  // Find source node
  R3SurfelNode *source_node = scene->Tree()->FindNodeByName(source_node_name);
  if (!source_node) {
    RNFail("Unable to find node with name %s\n", source_node_name);
    return 0;
  }

  // Create cluster objects 
  RNArray<R3SurfelObject *> *objects = CreateClusterObjects(scene, 
    source_node, NULL, parent_object, parent_node,                                                   
    max_neighbors, max_neighbor_distance, 
    max_offplane_distance, max_normal_angle, 
    min_points_per_object, chunk_size);
  if (!objects) {
    RNFail("No cluster objects created\n");
    return 0;
  }

  // Print statistics
  if (print_verbose) {
    printf("  Time = %.2f seconds\n", start_time.Elapsed());
    printf("  # Objects = %d\n", objects->NEntries());
    fflush(stdout);
  }

  // Delete array of objects
  delete objects;

  // Return success
  return 1;
}



static int
CreatePlanarObjects(R3SurfelScene *scene, 
  const char *parent_object_name, const char *parent_node_name, const char *source_node_name, 
  int max_neighbors, RNLength max_neighbor_distance, 
  RNLength max_offplane_distance, RNAngle max_normal_angle,
  RNArea min_area, RNScalar min_density, int min_points,
  RNLength grid_spacing, RNScalar accuracy_factor, RNLength chunk_size)
{
  // Start statistics
  RNTime start_time;
  start_time.Read();
  if (print_verbose) {
    printf("Creating planar objects ...\n");
    fflush(stdout);
  }

  // Find parent object
  R3SurfelObject *parent_object = scene->FindObjectByName(parent_object_name);
  if (!parent_object) {
    RNFail("Unable to find object with name %s\n", parent_object_name);
    return 0;
  }

  // Find parent node
  R3SurfelNode *parent_node = scene->Tree()->FindNodeByName(parent_node_name);
  if (!parent_node) {
    RNFail("Unable to find node with name %s\n", parent_node_name);
    return 0;
  }

  // Find source node
  R3SurfelNode *source_node = scene->Tree()->FindNodeByName(source_node_name);
  if (!source_node) {
    RNFail("Unable to find node with name %s\n", source_node_name);
    return 0;
  }

  // Create planar objects 
  RNArray<R3SurfelObject *> *objects = CreatePlanarObjects(scene, 
    source_node, NULL, parent_object, parent_node,                                                   
    max_neighbors, max_neighbor_distance, 
    max_offplane_distance, max_normal_angle,
    min_area, min_density, min_points, 
    grid_spacing, accuracy_factor, chunk_size);
  if (!objects) {
    RNFail("No planar objects created\n");
    return 0;
  }

  // Print statistics
  if (print_verbose) {
    printf("  Time = %.2f seconds\n", start_time.Elapsed());
    printf("  # Objects = %d\n", objects->NEntries());
    fflush(stdout);
  }

  // Delete array of objects
  delete objects;

  // Return success
  return 1;
}



////////////////////////////////////////////////////////////////////////
// OUTPUT FUNCTIONS
////////////////////////////////////////////////////////////////////////

static int
OutputBlobs(R3SurfelScene *scene, const char *directory_name)
{
  // Start statistics
  RNTime start_time;
  start_time.Read();
  if (print_verbose) {
    printf("Outputing blobs to %s ...\n", directory_name);
    fflush(stdout);
  }

  // Create directory
  char cmd[1024];
  sprintf(cmd, "mkdir -p %s", directory_name);
  system(cmd);

  // Output blob file for each object
  for (int i = 0; i < scene->NObjects(); i++) {
    R3SurfelObject *object = scene->Object(i);

    // Get object label
    R3SurfelLabel *label = object->GroundTruthLabel();
    if (!label) label = object->HumanLabel();
    int label_identifier = (label) ? label->Identifier() : 0;

    // Create pointset
    R3SurfelPointSet *pointset = object->PointSet();
    if (!pointset) continue;
    if (pointset->NPoints() == 0) { 
      delete pointset; 
      continue; 
    }

    // Create filename
    char filename[1024];
    R3Point centroid = object->Centroid();
    sprintf(filename, "%s/%d_%.3f_%.3f_%.3f.xyz", directory_name,
      label_identifier, centroid.X(), centroid.Y(), centroid.Z());

    // Open file
    FILE *fp = fopen(filename, "w");
    if (!fp) {
      RNFail("Unable to open xyz file %s\n", filename);
      delete pointset;
      return 0;
    }

    // Write points to file
    for (int j = 0; j < pointset->NPoints(); j++) {
      R3SurfelPoint *point = pointset->Point(j);
      R3Point position = point->Position();
      // RNRgb rgb = point->Rgb();
      fprintf(fp, "%.6f %.6f %.6f\n", position.X(), position.Y(), position.Z());
    }

    // Close file
    fclose(fp);

    // Print message
    if (print_debug) {
      printf("%3d %8.3f %8.3f %8.3f : %6d %g\n", label_identifier, 
        centroid.X(), centroid.Y(), centroid.Z(),
        pointset->NPoints(), object->Complexity());
    }

    // Delete pointset
    delete pointset;
  }

  // Print statistics
  if (print_verbose) {
    printf("  Time = %.2f seconds\n", start_time.Elapsed());
    printf("  # Blobs = %d\n", scene->NObjects());
    fflush(stdout);
  }
  
  // Return success
  return 1;
}



////////////////////////////////////////////////////////////////////////
// PROGRAM ARGUMENT PARSING
////////////////////////////////////////////////////////////////////////

static int
CheckForNumber(const char *str)
{
  // Check if string is a number
  for (const char *strp = str; *strp; strp++) {
    if (isdigit(*strp)) continue;
    if (*strp == '-') continue;
    if (*strp == '-') continue;
    if (*strp == '+') continue;
    if (*strp == '.') continue;
    return 0;
  }

  // Passed all tests
  return 1;
}



static int 
CheckForArgument(int argc, char **argv, const char *argument)
{
  // Check for -v in program arguments
  for (int i = 1; i < argc; i++) 
    if (!strcmp(argv[i], argument)) return 1;
  return 0;
}



R3SurfelConstraint *
ParseConstraint(int& argc, char **& argv)
{
  // Check arguments
  if (argc == 0) return NULL;

  // Check constraint type
  argc--; argv++;
  const char *constraint_type = *argv;
  if (!strcmp(constraint_type, "BoundingBox")) {
    // Read bounding box coordinates
    argc--; argv++; double x1 = atof(*argv);
    argc--; argv++; double y1 = atof(*argv);
    argc--; argv++; double z1 = atof(*argv);
    argc--; argv++; double x2 = atof(*argv);
    argc--; argv++; double y2 = atof(*argv);
    argc--; argv++; double z2 = atof(*argv);
    R3Box box(x1, y1, z1, x2, y2, z2);
    R3SurfelBoxConstraint *constraint = new R3SurfelBoxConstraint(box);
    return constraint;
  }
  else if (!strcmp(constraint_type, "OverheadGrid")) {
    // Read overhead grid
    argc--; argv++; const char *overhead_grid_name = *argv;
    R2Grid *overhead_grid = ReadGrid(overhead_grid_name);
    if (!overhead_grid) return NULL;

    // Parse comparison type
    argc--; argv++; const char *comparison_type_string = *argv;
    int comparison_type = 0;
    if (!strcmp(comparison_type_string, "NotEqual")) comparison_type = R3_SURFEL_CONSTRAINT_NOT_EQUAL;
    else if (!strcmp(comparison_type_string, "Equal")) comparison_type = R3_SURFEL_CONSTRAINT_EQUAL;
    else if (!strcmp(comparison_type_string, "Greater")) comparison_type = R3_SURFEL_CONSTRAINT_GREATER;
    else if (!strcmp(comparison_type_string, "GreaterOrEqual")) comparison_type = R3_SURFEL_CONSTRAINT_GREATER_OR_EQUAL;
    else if (!strcmp(comparison_type_string, "Less")) comparison_type = R3_SURFEL_CONSTRAINT_LESS;
    else if (!strcmp(comparison_type_string, "LessOrEqual")) comparison_type = R3_SURFEL_CONSTRAINT_LESS_OR_EQUAL;
    else { RNFail("Unrecognized constraint comparison type: %s\n", comparison_type_string); return NULL; }

    // Parse surfel operand
    argc--; argv++; const char *surfel_operand_string = *argv;
    RNScalar surfel_operand_value = 0;
    int surfel_operand_type = R3_SURFEL_CONSTRAINT_OPERAND;
    if (!strcmp(surfel_operand_string, "X")) surfel_operand_type = R3_SURFEL_CONSTRAINT_X;
    else if (!strcmp(surfel_operand_string, "Y")) surfel_operand_type = R3_SURFEL_CONSTRAINT_Y;
    else if (!strcmp(surfel_operand_string, "Z")) surfel_operand_type = R3_SURFEL_CONSTRAINT_Z;
    else if (CheckForNumber(surfel_operand_string)) surfel_operand_value = atof(surfel_operand_string); 
    else { RNFail("Unrecognized surfel operand: %s\n", surfel_operand_string); return NULL; }
  
    // Parse grid operand
    argc--; argv++; const char *grid_operand_string = *argv;
    RNScalar grid_operand_value = 0;
    int grid_operand_type = R3_SURFEL_CONSTRAINT_OPERAND;
    if (!strcmp(grid_operand_string, "Value")) grid_operand_type = R3_SURFEL_CONSTRAINT_VALUE;
    else if (CheckForNumber(grid_operand_string)) grid_operand_value = atof(grid_operand_string); 
    else { RNFail("Unrecognized grid operand: %s\n", grid_operand_string); return NULL; }

    // Parse epsilon
    argc--; argv++; double epsilon = atof(*argv);

    // Create constraint
    R3SurfelOverheadGridConstraint *constraint = new R3SurfelOverheadGridConstraint(
      overhead_grid, comparison_type, 
      surfel_operand_type, grid_operand_type, 
      surfel_operand_value, grid_operand_value, 
      epsilon);
    return constraint;
  }

  // Did not recognize constraint type
  RNFail("Unrecognized constraint type: %s\n", constraint_type); 
  return NULL;
}



////////////////////////////////////////////////////////////////////////
// Main
////////////////////////////////////////////////////////////////////////

int main(int argc, char **argv)
{
  // Check program arguments
  if (argc < 3) {
    // Print usage
    RNFail("Usage: surfelprocess scenefile databasefile [operations]\n");
    exit(-1);
  }

  // Parse program arguments
  scene_name = argv[1];
  database_name = argv[2];
  aerial_only = CheckForArgument(argc, argv, "-aerial_only");
  terrestrial_only = CheckForArgument(argc, argv, "-terrestrial_only");
  print_verbose = CheckForArgument(argc, argv, "-v");
  print_debug = CheckForArgument(argc, argv, "-debug");

  // Open scene
  R3SurfelScene *scene = OpenScene(scene_name, database_name);
  if (!scene) exit(-1);

  // Start statistics
  RNTime start_time;
  start_time.Read();
  if (print_verbose) {
    printf("Processing scene ...\n");
    fflush(stdout);
  }
  
  // Execute operations
  int noperations = 0;
  argc -= 3; argv += 3;
  while (argc > 0) {
    if (!strcmp(*argv, "-v")) print_verbose = 1;
    else if (!strcmp(*argv, "-debug")) print_debug = 1;
    else if (!strcmp(*argv, "-aerial_only")) aerial_only = 1;
    else if (!strcmp(*argv, "-terrestrial_only")) terrestrial_only = 1;
    else if (!strcmp(*argv, "-create_node")) { 
      argc--; argv++; char *node_name = *argv; 
      argc--; argv++; char *parent_name = *argv; 
      if (!CreateNode(scene, node_name, parent_name)) exit(-1);
      noperations++;
    }
    else if (!strcmp(*argv, "-create_object")) { 
      argc--; argv++; char *object_name = *argv; 
      argc--; argv++; char *parent_name = *argv; 
      argc--; argv++; char *node_name = *argv; 
      if (!CreateObject(scene, object_name, parent_name, node_name)) exit(-1);
      noperations++;
    }
    else if (!strcmp(*argv, "-create_label")) { 
      argc--; argv++; char *label_name = *argv; 
      argc--; argv++; char *parent_name = *argv; 
      if (!CreateLabel(scene, label_name, parent_name)) exit(-1);
      noperations++;
    }
    else if (!strcmp(*argv, "-load_surfels")) { 
      argc--; argv++; char *surfels_filename = *argv; 
      argc--; argv++; char *node_name = *argv; 
      argc--; argv++; char *parent_node_name = *argv; 
      if (!LoadSurfels(scene, surfels_filename, NULL, NULL, 
        node_name, parent_node_name)) exit(-1);
      noperations++;
    }
    else if (!strcmp(*argv, "-load_surfels_list")) { 
      argc--; argv++; char *list_filename = *argv; 
      argc--; argv++; char *parent_node_name = *argv; 
      if (!LoadSurfelsList(scene, list_filename, 
        NULL, parent_node_name)) exit(-1);
      noperations++;
    }
    else if (!strcmp(*argv, "-load_object")) { 
      argc--; argv++; char *surfels_filename = *argv; 
      argc--; argv++; char *object_name = *argv; 
      argc--; argv++; char *parent_object_name = *argv; 
      argc--; argv++; char *node_name = *argv; 
      argc--; argv++; char *parent_node_name = *argv; 
      if (!LoadSurfels(scene, surfels_filename, 
        object_name, parent_object_name, 
        node_name, parent_node_name)) exit(-1);
      noperations++;
    }
    else if (!strcmp(*argv, "-load_object_list")) { 
      argc--; argv++; char *list_filename = *argv; 
      argc--; argv++; char *parent_object_name = *argv; 
      argc--; argv++; char *parent_node_name = *argv; 
      if (!LoadSurfelsList(scene, list_filename, 
        parent_object_name, parent_node_name)) exit(-1);
      noperations++;
    }
    else if (!strcmp(*argv, "-load_label_list")) { 
      argc--; argv++; char *list_filename = *argv; 
      argc--; argv++; char *parent_label_name = *argv; 
      if (!LoadLabelList(scene, list_filename, parent_label_name)) exit(-1);
      noperations++;
    }
    else if (!strcmp(*argv, "-load_assignment_list")) { 
      argc--; argv++; char *list_filename = *argv; 
      if (!LoadAssignmentList(scene, list_filename)) exit(-1);
      noperations++;
    }
    else if (!strcmp(*argv, "-load_feature_list")) { 
      argc--; argv++; char *list_filename = *argv; 
      if (!LoadFeatureList(scene, list_filename)) exit(-1);
      noperations++;
    }
    else if (!strcmp(*argv, "-load_mesh")) { 
      argc--; argv++; char *mesh_filename = *argv; 
      argc--; argv++; char *parent_object_name = *argv; 
      argc--; argv++; char *parent_node_name = *argv; 
      argc--; argv++; double surfel_spacing = atof(*argv); 
      if (!LoadSurfelsFromMesh(scene, mesh_filename,
        parent_object_name, parent_node_name,
        surfel_spacing)) exit(-1);
      noperations++;
    }
    else if (!strcmp(*argv, "-load_scene")) { 
      argc--; argv++; char *scene_filename = *argv; 
      argc--; argv++; char *database_filename = *argv; 
      argc--; argv++; char *parent_object_name = *argv; 
      argc--; argv++; char *parent_label_name = *argv; 
      argc--; argv++; char *parent_node_name = *argv; 
      if (!LoadScene(scene, scene_filename, database_filename, 
        parent_object_name, parent_label_name, parent_node_name)) exit(-1);
      noperations++;
    }
    // else if (!strcmp(*argv, "-transfer_labels")) { 
    //   argc--; argv++; char *label_scene_filename = *argv; 
    //   argc--; argv++; char *label_database_filename = *argv; 
    //   if (!TransferLabels(scene, label_scene_filename, label_database_filename)) exit(-1);
    //     noperations++;
    // }
    else if (!strcmp(*argv, "-mask")) { 
      argc--; argv++; char *source_node_name = *argv; 
      R3SurfelConstraint *constraint = ParseConstraint(argc, argv);
      if (!constraint) exit(-1);
      if (!Mask(scene, source_node_name, constraint)) exit(-1);
      delete constraint;
      noperations++;
    }
    else if (!strcmp(*argv, "-remove_objects")) { 
      if (!RemoveObjects(scene)) exit(-1);
      noperations++;
    }
    else if (!strcmp(*argv, "-remove_labels")) { 
      if (!RemoveLabels(scene)) exit(-1);
      noperations++;
    }
    else if (!strcmp(*argv, "-remove_interior_nodes")) { 
      if (!RemoveInteriorNodes(scene)) exit(-1);
      noperations++;
    }
    else if (!strcmp(*argv, "-estimate_surfel_colors")) { 
      argc--; argv++; char *image_directory = *argv; 
      argc--; argv++; double depth_scale = atof(*argv); 
      argc--; argv++; double depth_exponent = atof(*argv); 
      if (!ReadImageDirectory(scene, image_directory, depth_scale, depth_exponent)) exit(-1);
      if (!EstimateSurfelColors(scene)) exit(-1);
      noperations++;
    }
    else if (!strcmp(*argv, "-order_surfel_identifiers")) { 
      if (!OrderSurfelIdentifiers(scene)) exit(-1);
      noperations++;
    }
    else if (!strcmp(*argv, "-transform_with_configuration_file")) { 
      argc--; argv++; char *configuration_filename = *argv; 
      if (!TransformWithConfigurationFile(scene, configuration_filename)) exit(-1);
      noperations++;
    }
    else if (!strcmp(*argv, "-inverse_transform_with_configuration_file")) { 
      argc--; argv++; char *configuration_filename = *argv; 
      if (!TransformWithConfigurationFile(scene, configuration_filename, TRUE)) exit(-1);
      noperations++;
    }
    else if (!strcmp(*argv, "-create_object_relationships")) { 
      argc--; argv++; double max_gap_distance = atof(*argv);
      argc--; argv++; double max_plane_offset = atof(*argv);
      argc--; argv++; double max_normal_angle = atof(*argv);
      argc--; argv++; double min_overlap = atof(*argv);
      if (!CreateObjectRelationships(scene, max_gap_distance, max_plane_offset, max_normal_angle, min_overlap)) exit(-1);
      noperations++;
    }
    else if (!strcmp(*argv, "-create_cluster_objects")) { 
      argc--; argv++; char *parent_object_name = *argv; 
      argc--; argv++; char *parent_node_name = *argv; 
      argc--; argv++; char *source_node_name = *argv; 
      argc--; argv++; int max_neighbors = atoi(*argv); 
      argc--; argv++; double max_neighbor_distance = atof(*argv); 
      argc--; argv++; double max_offplane_distance = atof(*argv); 
      argc--; argv++; double max_normal_angle = atof(*argv); 
      argc--; argv++; int min_points_per_object = atoi(*argv); 
      argc--; argv++; double chunk_size= atof(*argv); 
      if (!CreateClusterObjects(scene, parent_object_name, parent_node_name, source_node_name,
        chunk_size, max_neighbors, max_neighbor_distance, 
        max_offplane_distance, max_normal_angle, min_points_per_object)) exit(-1);
      noperations++;
    }
    else if (!strcmp(*argv, "-create_planar_objects")) { 
      argc--; argv++; char *parent_object_name = *argv; 
      argc--; argv++; char *parent_node_name = *argv; 
      argc--; argv++; char *source_node_name = *argv; 
      argc--; argv++; int max_neighbors = atoi(*argv); 
      argc--; argv++; double max_neighbor_distance = atof(*argv); 
      argc--; argv++; double max_offplane_distance = atof(*argv); 
      argc--; argv++; double max_normal_angle = atof(*argv); 
      argc--; argv++; double min_area = atof(*argv); 
      argc--; argv++; double min_density = atof(*argv); 
      argc--; argv++; double min_points = atof(*argv); 
      argc--; argv++; double grid_spacing = atof(*argv); 
      argc--; argv++; double accuracy_factor = atof(*argv); 
      argc--; argv++; double chunk_size= atof(*argv); 
      if (!CreatePlanarObjects(scene, parent_object_name, parent_node_name, source_node_name,
        chunk_size, max_neighbors, max_neighbor_distance, max_offplane_distance, max_normal_angle,
        min_area, min_density, min_points, grid_spacing, accuracy_factor)) exit(-1);
      noperations++;
    }
    else if (!strcmp(*argv, "-create_multiresolution_hierarchy")) { 
      const char *node_name = "Root";
      int max_parts_per_node = 8;
      int max_blocks_per_node = 32;
      RNScalar max_node_complexity = 1024;
      RNScalar max_block_complexity = 1024;
      RNLength max_leaf_extent = 10.0;
      RNLength max_block_extent = 10.0;
      RNScalar multiresolution_factor = 0.25;
      int max_levels = 64;
      if (!SplitSurfelTreeNodes(scene, node_name, 
        max_parts_per_node, max_blocks_per_node, 
        max_node_complexity, max_block_complexity,
        max_leaf_extent, max_block_extent, max_levels)) exit(-1);
      if (!CreateMultiresolutionBlocks(scene, node_name, 
        multiresolution_factor, max_node_complexity)) exit(-1);
      noperations++;

    }
    else if (!strcmp(*argv, "-create_tree_hierarchy")) { 
      argc--; argv++; char *node_name = *argv; 
      argc--; argv++; int max_parts_per_node = atoi(*argv);
      argc--; argv++; int max_blocks_per_node = atoi(*argv);
      argc--; argv++; double max_node_complexity = atof(*argv); 
      argc--; argv++; double max_block_complexity = atof(*argv); 
      argc--; argv++; double max_leaf_extent = atof(*argv); 
      argc--; argv++; double max_block_extent = atof(*argv); 
      argc--; argv++; double multiresolution_factor = atof(*argv); 
      argc--; argv++; int max_levels = atoi(*argv); 
      if (!SplitSurfelTreeNodes(scene, node_name, 
        max_parts_per_node, max_blocks_per_node, 
        max_node_complexity, max_block_complexity,
        max_leaf_extent, max_block_extent, max_levels)) exit(-1);
      if (!CreateMultiresolutionBlocks(scene, node_name, 
        multiresolution_factor, max_node_complexity)) exit(-1);
      noperations++;
    }
    else if (!strcmp(*argv, "-split_nodes")) { 
      argc--; argv++; char *node_name = *argv; 
      argc--; argv++; int max_parts_per_node = atoi(*argv);
      argc--; argv++; int max_blocks_per_node = atoi(*argv);
      argc--; argv++; double max_node_complexity = atof(*argv); 
      argc--; argv++; double max_block_complexity = atof(*argv); 
      argc--; argv++; double max_leaf_extent = atof(*argv); 
      argc--; argv++; double max_block_extent = atof(*argv);
      argc--; argv++; int max_levels = atoi(*argv); 
      if (!SplitSurfelTreeNodes(scene, node_name, 
        max_parts_per_node, max_blocks_per_node, 
        max_node_complexity, max_block_complexity,
        max_leaf_extent, max_block_extent, max_levels)) exit(-1);
      noperations++;
    }
    else if (!strcmp(*argv, "-create_multiresolution_nodes")) { 
      argc--; argv++; char *node_name = *argv; 
      argc--; argv++; double min_complexity = atof(*argv); 
      argc--; argv++; double min_resolution = atof(*argv); 
      argc--; argv++; double min_multiresolution_factor = atof(*argv); 
      if (!CreateMultiresolutionNodes(scene, node_name, 
        min_complexity, min_resolution, min_multiresolution_factor)) exit(-1);
      noperations++;
    }
    else if (!strcmp(*argv, "-create_multiresolution_blocks")) { 
      argc--; argv++; char *node_name = *argv; 
      argc--; argv++; double multiresolution_factor = atof(*argv); 
      argc--; argv++; double max_node_complexity = atof(*argv); 
      if (!CreateMultiresolutionBlocks(scene, node_name, 
        multiresolution_factor, max_node_complexity)) exit(-1);
      noperations++;
    }
    else if (!strcmp(*argv, "-output_blobs")) { 
      argc--; argv++; char *blob_directory_name = *argv; 
      if (!OutputBlobs(scene, blob_directory_name)) exit(-1);
      noperations++;
    }
    else { 
      RNFail("Invalid operation: %s", *argv); 
      exit(1); 
    }
    argv++; argc--;
  }

  // Print statistics
  if (print_verbose) {
    printf("  Time = %.2f\n", start_time.Elapsed());
    printf("  # Operations = %d\n", noperations);
    fflush(stdout);
  }
  
  // Close scene
  if (!CloseScene(scene)) exit(-1);

  // Delete scene
  delete scene;
  
  // Return success 
  return 0;
}



