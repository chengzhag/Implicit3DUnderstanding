// Source file for the scene converter program



// Include files 

namespace gaps {}
using namespace gaps;
#include "R3Graphics/R3Graphics.h"



// Program arguments

static const char *input_name = NULL;
static const char *output_name = NULL;
static char *input_categories_name = NULL;
static char *input_lights_name = NULL;
static const char *select_nodes_in_subtree_name = NULL;
static const char *select_nodes_in_grid_filename = NULL;
static R3Box select_nodes_in_bbox(FLT_MAX, FLT_MAX, FLT_MAX, -FLT_MAX, -FLT_MAX, -FLT_MAX);
static R3Affine xform(R4Matrix(1,0,0,0, 0,1,0,0, 0,0,1,0, 0,0,0,1), 0);
static int remove_references = 0;
static int remove_hierarchy = 0;
static int remove_transformations = 0;
static RNLength max_edge_length = 0;
static int print_verbose = 0;



////////////////////////////////////////////////////////////////////////
// I/O STUFF
////////////////////////////////////////////////////////////////////////

static R3Scene *
ReadScene(const char *filename)
{
  // Start statistics
  RNTime start_time;
  start_time.Read();

  // Allocate scene
  R3Scene *scene = new R3Scene();
  if (!scene) {
    RNFail("Unable to allocate scene for %s\n", filename);
    return NULL;
  }

  // Read scene from file
  if (!scene->ReadFile(filename)) {
    delete scene;
    return NULL;
  }

  // Print statistics
  if (print_verbose) {
    printf("Read scene from %s ...\n", filename);
    printf("  Time = %.2f seconds\n", start_time.Elapsed());
    printf("  # Nodes = %d\n", scene->NNodes());
    printf("  # Lights = %d\n", scene->NLights());
    printf("  # Materials = %d\n", scene->NMaterials());
    printf("  # Brdfs = %d\n", scene->NBrdfs());
    printf("  # Textures = %d\n", scene->NTextures());
    printf("  # Referenced models = %d\n", scene->NReferencedScenes());
    fflush(stdout);
  }

  // Return scene
  return scene;
}



static int
WriteScene(R3Scene *scene, const char *filename)
{
  // Start statistics
  RNTime start_time;
  start_time.Read();

  // Write scene from file
  if (!scene->WriteFile(filename)) return 1;

  // Print statistics
  if (print_verbose) {
    printf("Wrote scene to %s ...\n", filename);
    printf("  Time = %.2f seconds\n", start_time.Elapsed());
    printf("  # Nodes = %d\n", scene->NNodes());
    printf("  # Lights = %d\n", scene->NLights());
    printf("  # Materials = %d\n", scene->NMaterials());
    printf("  # Brdfs = %d\n", scene->NBrdfs());
    printf("  # Textures = %d\n", scene->NTextures());
    printf("  # Referenced models = %d\n", scene->NReferencedScenes());
    fflush(stdout);
  }

  // Return success
  return 1;
}



static R3Grid *
ReadGrid(const char *input_name)
{
  // Start statistics
  RNTime start_time;
  start_time.Read();

  // Allocate a grid
  R3Grid *grid = new R3Grid();
  if (!grid) {
    RNFail("Unable to allocate grid");
    return NULL;
  }

  // Read grid
  int status = grid->ReadFile(input_name);
  if (!status) {
    RNFail("Unable to read grid file %s", input_name);
    return NULL;
  }

  // Print statistics
  if (print_verbose) {
    printf("Read grid ...\n");
    printf("  Time = %.2f seconds\n", start_time.Elapsed());
    printf("  Resolution = %d %d %d\n", grid->XResolution(), grid->YResolution(), grid->ZResolution());
    printf("  Spacing = %g\n", grid->GridToWorldScaleFactor());
    printf("  Cardinality = %d\n", grid->Cardinality());
    printf("  Volume = %g\n", grid->Volume());
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



static int
ReadCategories(R3Scene *scene, const char *filename)
{
  // Start statistics
  RNTime start_time;
  start_time.Read();

  // Read file
  if (!scene->ReadSUNCGModelFile(filename)) return 0;

  // Print statistics
  if (print_verbose) {
    printf("Read categories from %s ...\n", filename);
    printf("  Time = %.2f seconds\n", start_time.Elapsed());
    fflush(stdout);
  }

  // Return success
  return 1;
} 



static int
ReadLights(R3Scene *scene, const char *filename)
{
  // Start statistics
  RNTime start_time;
  start_time.Read();

  // Read lights file
  if (!scene->ReadSUNCGLightsFile(filename)) return 0;

  // Print statistics
  if (print_verbose) {
    printf("Read lights from %s ...\n", filename);
    printf("  Time = %.2f seconds\n", start_time.Elapsed());
    printf("  # Lights = %d\n", scene->NLights());
    fflush(stdout);
  }

  // Return success
  return 1;
}
  


static int
ReadMatrix(R4Matrix& m, const char *filename)
{
  // Open file
  FILE *fp = fopen(filename, "r");
  if (!fp) {
    RNFail("Unable to open matrix file: %s\n", filename);
    return 0;
  }

  // Read matrix from file
  for (int i = 0; i < 4; i++) {
    for (int j = 0; j < 4; j++) {
      double value;
      fscanf(fp, "%lf", &value);
      m[i][j] = value;
    }
  }

  // Close file
  fclose(fp);

  // Return success
  return 1;
}



////////////////////////////////////////////////////////////////////////
// PROCESSING
////////////////////////////////////////////////////////////////////////

static void
CountOverlapAreas(R3SceneNode *node, R3Grid *grid, const R3Affine& parent_transformation, RNArea& overlap_area, RNArea& total_area)
{
  // Update transformation
  R3Affine transformation = R3identity_affine;
  transformation.Transform(parent_transformation);
  transformation.Transform(node->Transformation());

  // Count overlapping triangle areas in references
  for (int i = 0; i < node->NReferences(); i++) {
    R3SceneReference *reference = node->Reference(i);
    R3Scene *referenced_scene = reference->ReferencedScene();
    CountOverlapAreas(referenced_scene->Root(), grid, transformation, overlap_area, total_area);
  }

  // Count overlapping triangle areas in elements
  for (int i = 0; i < node->NElements(); i++) {
    R3SceneElement *element = node->Element(i);
    for (int j = 0; j < element->NShapes(); j++) {
      R3Shape *shape = element->Shape(j);
      if (shape->ClassID() == R3TriangleArray::CLASS_ID()) {
        R3TriangleArray *triangles = (R3TriangleArray *) shape;
        for (int k = 0; k < triangles->NTriangles(); k++) {
          R3Triangle *triangle = triangles->Triangle(k);
          R3Point centroid = triangle->Centroid();
          centroid.Transform(transformation);
          RNScalar grid_value = grid->WorldValue(centroid);
          RNArea area = triangle->Area();
          if (grid_value > RN_EPSILON) overlap_area += area;
          total_area += area;
        }
      }
    }
  }
  
  // Count overlapping triangle areas in children
  for (int i = 0; i < node->NChildren(); i++) {
    R3SceneNode *child = node->Child(i);
    CountOverlapAreas(child, grid, transformation, overlap_area, total_area);
  }
}



static RNScalar
Overlap(R3SceneNode *node, R3Grid *grid)
{
  RNScalar overlap_area = 0;
  RNScalar total_area = 0;
  R3Affine parent_transformation = node->CumulativeParentTransformation();
  CountOverlapAreas(node, grid, parent_transformation, overlap_area, total_area);
  if (total_area < RN_EPSILON) return 0;
  return overlap_area / total_area;
}



static void
RemoveNodes(R3Scene *scene, R3SceneNode *node,
  R3SceneNode *select_subtree_node, R3Grid *select_grid, const R3Box& select_bbox)
{
  // Copy array of children -- because will be edited as iterate
  RNArray<R3SceneNode *> children;
  for (int i = 0; i < node->NChildren(); i++) {
    R3SceneNode *child = node->Child(i);
    children.Insert(child);
  }

  // Visit children recursively 
  for (int i = 0; i < children.NEntries(); i++) {
    R3SceneNode *child = children.Kth(i);
    RemoveNodes(scene, child, select_subtree_node, select_grid, select_bbox);
  }

  // Check if node must remain
  if (node == scene->Root()) return;
  if (node->NChildren() > 0) return;

  // Innocent until proven guilty
  RNBoolean remove = FALSE;

  // Check if should remove based on emptiness
  if ((node->NChildren() == 0) && (node->NReferences() == 0) && (node->NElements() == 0)) {
    remove = TRUE;
  }

  // Check if should remove based on bbox selection
  if (!select_bbox.IsEmpty()) {
    if (!R3Intersects(node->WorldBBox(), select_bbox)) {
      remove = TRUE;
    }
  }

  // Check if should remove based on subtree selection
  if (select_subtree_node) {
    if (node != select_subtree_node) {
      if (!node->IsAncestor(select_subtree_node)) { 
        if (!node->IsDecendent(select_subtree_node)) {
          remove = TRUE;
        }
      }
    }
  }

  // Check if should remove based on grid selection
  if (select_grid) {
    RNArea min_overlap = 0.5;
    RNScalar overlap = Overlap(node, select_grid);
    if (overlap < min_overlap) remove = TRUE;
    // R3Point centroid = node->WorldBBox().Centroid();
    // RNScalar grid_value = select_grid->WorldValue(centroid);
    // if (RNIsNegativeOrZero(grid_value)) remove = TRUE;
  }
  
  // Delete node 
  if (remove) delete node;
}



static int
RemoveNodes(R3Scene *scene)
{
  // Find select subtree node
  R3SceneNode *select_subtree_node = NULL;
  if (select_nodes_in_subtree_name) {
    select_subtree_node = scene->Node(select_nodes_in_subtree_name);
    if (!select_subtree_node) {
      RNFail("Unable to find select subtree node %s\n", select_nodes_in_subtree_name);
      return 0;
    }
  }

  // Read select grid
  R3Grid *select_grid = NULL;
  if (select_nodes_in_grid_filename) {
    select_grid = ReadGrid(select_nodes_in_grid_filename);
    if (!select_grid) return 0;
  }

  // Remove nodes recursively in post order
  RemoveNodes(scene, scene->Root(), select_subtree_node, select_grid, select_nodes_in_bbox);

  // Return success
  return 1;
}



////////////////////////////////////////////////////////////////////////
// PROGRAM ARGUMENT PARSING
////////////////////////////////////////////////////////////////////////

static int
ParseArgs(int argc, char **argv)
{
  // Parse arguments
  argc--; argv++;
  while (argc > 0) {
    if ((*argv)[0] == '-') {
      R3Affine prev_xform = xform;
      if (!strcmp(*argv, "-v")) print_verbose = 1;
      else if (!strcmp(*argv, "-remove_references")) remove_references = 1;
      else if (!strcmp(*argv, "-remove_hierarchy")) remove_hierarchy = 1;
      else if (!strcmp(*argv, "-remove_transformations")) remove_transformations = 1;
      else if (!strcmp(*argv, "-select_nodes_in_subtree")) {
        argv++; argc--; select_nodes_in_subtree_name = *argv;
      }
      else if (!strcmp(*argv, "-select_nodes_in_grid")) {
        argv++; argc--; select_nodes_in_grid_filename = *argv;
      }
      else if (!strcmp(*argv, "-select_nodes_in_bbox")) {
        argv++; argc--; select_nodes_in_bbox[0][0] = atof(*argv);
        argv++; argc--; select_nodes_in_bbox[0][1] = atof(*argv);
        argv++; argc--; select_nodes_in_bbox[0][2] = atof(*argv);
        argv++; argc--; select_nodes_in_bbox[1][0] = atof(*argv);
        argv++; argc--; select_nodes_in_bbox[1][1] = atof(*argv);
        argv++; argc--; select_nodes_in_bbox[1][2] = atof(*argv);
      }
      else if (!strcmp(*argv, "-scale")) {
        argv++; argc--; xform.Scale(atof(*argv));
      }
      else if (!strcmp(*argv, "-tx")) {
        argv++; argc--; xform = R3identity_affine; xform.XTranslate(atof(*argv)); xform.Transform(prev_xform);
      }
      else if (!strcmp(*argv, "-ty")) {
        argv++; argc--; xform = R3identity_affine; xform.YTranslate(atof(*argv)); xform.Transform(prev_xform);
      }
      else if (!strcmp(*argv, "-tz")) {
        argv++; argc--; xform = R3identity_affine; xform.ZTranslate(atof(*argv)); xform.Transform(prev_xform);
      }
      else if (!strcmp(*argv, "-sx")) {
        argv++; argc--; xform = R3identity_affine; xform.XScale(atof(*argv)); xform.Transform(prev_xform);
      }
      else if (!strcmp(*argv, "-sy")) {
        argv++; argc--; xform = R3identity_affine; xform.YScale(atof(*argv)); xform.Transform(prev_xform);
      }
      else if (!strcmp(*argv, "-sz")) {
        argv++; argc--; xform = R3identity_affine; xform.ZScale(atof(*argv)); xform.Transform(prev_xform);
      }
      else if (!strcmp(*argv, "-rx")) {
        argv++; argc--; xform = R3identity_affine; xform.XRotate(RN_PI*atof(*argv)/180.0); xform.Transform(prev_xform);
      }
      else if (!strcmp(*argv, "-ry")) {
        argv++; argc--; xform = R3identity_affine; xform.YRotate(RN_PI*atof(*argv)/180.0); xform.Transform(prev_xform);
      }
      else if (!strcmp(*argv, "-rz")) {
        argv++; argc--; xform = R3identity_affine; xform.ZRotate(RN_PI*atof(*argv)/180.0); xform.Transform(prev_xform);
      }
      else if (!strcmp(*argv, "-xform")) {
        argv++; argc--; R4Matrix m;  if (ReadMatrix(m, *argv)) { xform = R3identity_affine; xform.Transform(R3Affine(m)); xform.Transform(prev_xform);}
      } 
      else if (!strcmp(*argv, "-max_edge_length")) {
        argv++; argc--; max_edge_length = atof(*argv);
      }
      else if (!strcmp(*argv, "-categories")) {
        argc--; argv++; input_categories_name = *argv;
      }
      else if (!strcmp(*argv, "-lights")) {
        argv++; argc--; input_lights_name = *argv;
      }
      else {
        RNFail("Invalid program argument: %s", *argv);
        exit(1);
      }
      argv++; argc--;
    }
    else {
      if (!input_name) input_name = *argv;
      else if (!output_name) output_name = *argv;
      else { RNFail("Invalid program argument: %s", *argv); exit(1); }
      argv++; argc--;
    }
  }

  // Check input filename
  if (!input_name || !output_name) {
    RNFail("Usage: scn2scn inputfile outputfile [options]\n");
    return 0;
  }

  // Return OK status 
  return 1;
}



////////////////////////////////////////////////////////////////////////
// MAIN
////////////////////////////////////////////////////////////////////////

int main(int argc, char **argv)
{
  // Check number of arguments
  if (!ParseArgs(argc, argv)) exit(1);

  // Read scene
  R3Scene *scene = ReadScene(input_name);
  if (!scene) exit(-1);

  // Read categories
  if (input_categories_name) {
    if (!ReadCategories(scene, input_categories_name)) exit(-1);
  }

  // Read lights
  if (input_lights_name) {
    if (!ReadLights(scene, input_lights_name)) exit(-1);
  }

  // Remove unselected nodes
  if (select_nodes_in_subtree_name || select_nodes_in_grid_filename || !select_nodes_in_bbox.IsEmpty()) {
    if (!RemoveNodes(scene)) exit(-1);
  }

  // Transform scene
  if (!xform.IsIdentity()) {
    R3SceneNode *root = scene->Root();
    R3Affine tmp = R3identity_affine;
    tmp.Transform(xform);
    tmp.Transform(root->Transformation());
    root->SetTransformation(tmp);
  }

  // Apply processing operations
  if (remove_references) scene->RemoveReferences();
  if (remove_hierarchy) scene->RemoveHierarchy();
  if (remove_transformations) scene->RemoveTransformations();
  if (max_edge_length > 0) scene->SubdivideTriangles(max_edge_length);

  // Write scene
  if (!WriteScene(scene, output_name)) exit(-1);

  // Return success 
  return 0;
}

















