// Source file for the pdb viewer program



////////////////////////////////////////////////////////////////////////
// Include files 
////////////////////////////////////////////////////////////////////////

namespace gaps {}
using namespace gaps;
#include "R3Graphics/R3Graphics.h"
#include "R3PointSet.h"
#include "fglut/fglut.h"



////////////////////////////////////////////////////////////////////////
// Global variables
////////////////////////////////////////////////////////////////////////

// Program variables
static char *image_name = NULL;
static RNScalar background_color[3] = { 1, 1, 1 };
static int print_verbose = 0;

// GLUT variables 
static int GLUTwindow = 0;
static int GLUTwindow_height = 1024;
static int GLUTwindow_width = 1024;
static int GLUTmouse[2] = { 0, 0 };
static int GLUTbutton[3] = { 0, 0, 0 };
static int GLUTmodifiers = 0;

// Application variables
static R3Viewer *viewer = NULL;
static R3Point initial_camera_origin = R3Point(0.0, 0.0, 0.0);
static R3Vector initial_camera_towards = R3Vector(0.0, 0.0, -1.0);
static R3Vector initial_camera_up = R3Vector(0.0, 1.0, 0.0);
static RNBoolean initial_camera = FALSE;
static R3Box world_bbox(0,0,0,0,0,0);
static R3Point world_origin(0, 0, 0);
static RNArray<R3Scene *> scenes;
static RNArray<R3Mesh *> meshes;
static RNArray<R3PointSet *> pointsets;
static RNArray<R3Grid *> grids;

// Display variables
static int show_scene_faces = 1;
static int show_scene_edges = 0;
static int show_grid_faces = 0;
static int show_grid_edges = 1;
static int show_grid_threshold = 0;
static int show_grid_slices[3] = { 0, 0, 0 };
static int show_mesh_faces = 1;
static int show_mesh_edges = 0;
static int show_mesh_vertices = 0;
static int show_points = 1;
static int show_bbox = 0;
static int show_axes = 0;
static int show_text_labels = 0;
static int show_normals = 0;
static int show_categories = 0;
static int show_segments = 0;
static int show_colors = 0;

// Model selection
static int nmodels = 0;
static const int max_models = 64;
static int selected_model_index = 0;

// Grid display variables
static RNScalar grid_thresholds[max_models] = { 0 };
static RNInterval grid_intervals[max_models] = { RNInterval(0,0) };
static int grid_slice_coords[3] = { 0, 0, 0 };



////////////////////////////////////////////////////////////////////////
// I/O functions
////////////////////////////////////////////////////////////////////////

static R3Scene *
ReadScene(char *filename)
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
    printf("  # Referenced scenes = %d\n", scene->NReferencedScenes());
    fflush(stdout);
  }

  // Return scene
  return scene;
}



static R3Mesh *
ReadMesh(const char *filename)
{
  // Start statistics
  RNTime start_time;
  start_time.Read();

  // Allocate mesh
  R3Mesh *mesh = new R3Mesh();
  if (!mesh) {
    RNFail("Unable to allocate mesh for %s\n", filename);
    return NULL;
  }

  // Read mesh from file
  if (!mesh->ReadFile(filename)) {
    delete mesh;
    return NULL;
  }

  // Print statistics
  if (print_verbose) {
    printf("Read mesh from %s ...\n", filename);
    printf("  Time = %.2f seconds\n", start_time.Elapsed());
    printf("  # Faces = %d\n", mesh->NFaces());
    printf("  # Edges = %d\n", mesh->NEdges());
    printf("  # Vertices = %d\n", mesh->NVertices());
    fflush(stdout);
  }

  // Return success
  return mesh;
}



static R3PointSet *
ReadPointSet(const char *filename)
{
  // Start statistics
  RNTime start_time;
  start_time.Read();

  // Allocate pointset
  R3PointSet *pointset = new R3PointSet();
  if (!pointset) {
    RNFail("Unable to allocate pointset for %s\n", filename);
    return NULL;
  }

  // Read pointset from file
  if (!pointset->ReadFile(filename)) {
    delete pointset;
    return NULL;
  }

  // Print statistics
  if (print_verbose) {
    printf("Read pointset from %s ...\n", filename);
    printf("  Time = %.2f seconds\n", start_time.Elapsed());
    printf("  # Points = %d\n", pointset->NPoints());
    fflush(stdout);
  }

  // Return success
  return pointset;
}



static R3Grid *
ReadGrid(char *filename)
{
  // Start statistics
  RNTime start_time;
  start_time.Read();

  // Allocate grid
  R3Grid *grid = new R3Grid();
  if (!grid) {
    RNFail("Unable to allocated grid");
    return NULL;
  }

  // Read grid 
  if (!grid->ReadFile(filename)) {
    RNFail("Unable to read grid file %s", filename);
    return NULL;
  }

  // Print statistics
  if (print_verbose) {
    const R3Box bbox = grid->WorldBox();
    RNInterval grid_range = grid->Range();
    printf("Read grid from %s ...\n", filename);
    printf("  Time = %.2f seconds\n", start_time.Elapsed());
    printf("  Resolution = %d %d %d\n", grid->XResolution(), grid->YResolution(), grid->ZResolution());
    printf("  World Box = ( %g %g %g ) ( %g %g %g )\n", bbox[0][0], bbox[0][1], bbox[0][2], bbox[1][0], bbox[1][1], bbox[1][2]);
    printf("  Spacing = %g\n", grid->GridToWorldScaleFactor());
    printf("  Cardinality = %d\n", grid->Cardinality());
    printf("  Volume = %g\n", grid->Volume());
    printf("  Minimum = %g\n", grid_range.Min());
    printf("  Maximum = %g\n", grid_range.Max());
    printf("  L1Norm = %g\n", grid->L1Norm());
    fflush(stdout);
  }

  // Return success
  return grid;
}



////////////////////////////////////////////////////////////////////////
// Coloring functions
////////////////////////////////////////////////////////////////////////

#if 0
static void
LoadIndex(int index, int tag)
{
  // Set color to represent an integer (24 bits)
  int k = index + 1;
  unsigned char color[4];
  color[0] = (k >> 16) & 0xFF;
  color[1] = (k >>  8) & 0xFF;
  color[2] = (k      ) & 0xFF;
  color[3] = tag;
  glColor4ubv(color);
}
#endif



static void
LoadColor(int k)
{
  // Make array of colors
  const int ncolors = 72;
  const RNRgb colors[ncolors] = {
    RNRgb(0.5, 0.5, 0.5), RNRgb(1, 0, 0), RNRgb(0, 0, 1), 
    RNRgb(0, 1, 0), RNRgb(0, 1, 1), RNRgb(1, 0, 1), 
    RNRgb(1, 0.5, 0), RNRgb(0, 1, 0.5), RNRgb(0.5, 0, 1), 
    RNRgb(0.5, 1, 0), RNRgb(0, 0.5, 1), RNRgb(1, 0, 0.5), 
    RNRgb(0.5, 0, 0), RNRgb(0, 0.5, 0), RNRgb(0, 0, 0.5), 
    RNRgb(0.5, 0.5, 0), RNRgb(0, 0.5, 0.5), RNRgb(0.5, 0, 0.5),
    RNRgb(0.7, 0, 0), RNRgb(0, 0.7, 0), RNRgb(0, 0, 0.7), 
    RNRgb(0.7, 0.7, 0), RNRgb(0, 0.7, 0.7), RNRgb(0.7, 0, 0.7), 
    RNRgb(0.7, 0.3, 0), RNRgb(0, 0.7, 0.3), RNRgb(0.3, 0, 0.7), 
    RNRgb(0.3, 0.7, 0), RNRgb(0, 0.3, 0.7), RNRgb(0.7, 0, 0.3), 
    RNRgb(0.3, 0, 0), RNRgb(0, 0.3, 0), RNRgb(0, 0, 0.3), 
    RNRgb(0.3, 0.3, 0), RNRgb(0, 0.3, 0.3), RNRgb(0.3, 0, 0.3),
    RNRgb(1, 0.3, 0.3), RNRgb(0.3, 1, 0.3), RNRgb(0.3, 0.3, 1), 
    RNRgb(1, 1, 0.3), RNRgb(0.3, 1, 1), RNRgb(1, 0.3, 1), 
    RNRgb(1, 0.5, 0.3), RNRgb(0.3, 1, 0.5), RNRgb(0.5, 0.3, 1), 
    RNRgb(0.5, 1, 0.3), RNRgb(0.3, 0.5, 1), RNRgb(1, 0.3, 0.5), 
    RNRgb(0.5, 0.3, 0.3), RNRgb(0.3, 0.5, 0.3), RNRgb(0.3, 0.3, 0.5), 
    RNRgb(0.5, 0.5, 0.3), RNRgb(0.3, 0.5, 0.5), RNRgb(0.5, 0.3, 0.5),
    RNRgb(0.3, 0.5, 0.5), RNRgb(0.5, 0.3, 0.5), RNRgb(0.5, 0.5, 0.3), 
    RNRgb(0.3, 0.3, 0.5), RNRgb(0.5, 0.3, 0.3), RNRgb(0.3, 0.5, 0.3), 
    RNRgb(0.3, 0.8, 0.5), RNRgb(0.5, 0.3, 0.8), RNRgb(0.8, 0.5, 0.3), 
    RNRgb(0.8, 0.3, 0.5), RNRgb(0.5, 0.8, 0.3), RNRgb(0.3, 0.5, 0.8), 
    RNRgb(0.8, 0.5, 0.5), RNRgb(0.5, 0.8, 0.5), RNRgb(0.5, 0.5, 0.8), 
    RNRgb(0.8, 0.8, 0.5), RNRgb(0.5, 0.8, 0.8), RNRgb(0.8, 0.5, 0.8)
  };

  // Load color
  if (k == -1) glColor3d(0.8, 0.8, 0.8);
  else if (k == 0) RNLoadRgb(colors[0]);
  else RNLoadRgb(colors[1 + (k % (ncolors-1))]);
}



////////////////////////////////////////////////////////////////////////
// Drawing functions
////////////////////////////////////////////////////////////////////////

static void
DrawText(const R3Point& p, const char *s)
{
  // Draw text string s and position p
  glRasterPos3d(p[0], p[1], p[2]);
  while (*s) glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12, *(s++));
}




static void 
DrawText(const R2Point& p, const char *s)
{
  // Draw text string s and position p
  R3Ray ray = viewer->WorldRay((int) p[0], (int) p[1]);
  R3Point position = ray.Point(2 * viewer->Camera().Near());
  glRasterPos3d(position[0], position[1], position[2]);
  while (*s) glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12, *(s++));
}



static void
DrawScene(R3Scene *scene)
{
  // Draw faces
  if (show_scene_faces) {
    glEnable(GL_LIGHTING);
    scene->LoadLights(2);
    scene->Draw();
  }

  // Draw edges
  if (show_scene_edges) {
    glDisable(GL_LIGHTING);
    glColor3d(0.0, 1.0, 0.0);
    scene->Draw(R3_EDGES_DRAW_FLAG);
  }
}



static void
DrawMesh(R3Mesh *mesh)
{
  // Draw mesh faces
  if (show_mesh_faces) {
    if (show_segments) glDisable(GL_LIGHTING);
    else if (show_categories) glDisable(GL_LIGHTING);
    else if (show_colors) glDisable(GL_LIGHTING);
    else { glEnable(GL_LIGHTING); glColor3d(0.8, 0.8, 0.8); }
    glBegin(GL_TRIANGLES);
    for (int i = 0; i < mesh->NFaces(); i++) {
      R3MeshFace *face = mesh->Face(i);
      if (show_segments) LoadColor(1 + mesh->FaceSegment(face)); 
      else if (show_categories) LoadColor(1 + mesh->FaceCategory(face));
      if (!show_colors) R3LoadNormal(mesh->FaceNormal(face));
      for (int j = 0; j < 3; j++) {
        R3MeshVertex *vertex = mesh->VertexOnFace(face, j);
        if (show_colors) R3LoadRgb(mesh->VertexColor(vertex));
        R3LoadPoint(mesh->VertexPosition(vertex));
      }
    }
    glEnd();
  }

  // Draw mesh edges
  if (show_mesh_edges) {
    glDisable(GL_LIGHTING);
    glColor3f(1.0, 0.0, 0.0);
    mesh->DrawEdges();
  }

  // Draw mesh vertices
  if (show_mesh_vertices) {
    if (show_colors) glDisable(GL_LIGHTING);
    else { glEnable(GL_LIGHTING); glColor3d(0.0, 0.5, 1.0); }
    glBegin(GL_POINTS);
    for (int i = 0; i < mesh->NVertices(); i++) {
      R3MeshVertex *vertex = mesh->Vertex(i);
      if (show_colors) R3LoadRgb(mesh->VertexColor(vertex));
      else R3LoadNormal(mesh->VertexNormal(vertex));
      R3LoadPoint(mesh->VertexPosition(vertex));
    }
    glEnd();
  }
}



static void
DrawPointSet(R3PointSet *pointset)
{
  // Draw points
  if (show_points) {
    glPointSize(3);
    glDisable(GL_LIGHTING);    
    if (!show_colors) glColor3d(0.0, 1.0, 0.0);
    glBegin(GL_POINTS);
    for (int i = 0; i < pointset->NPoints(); i++) {
      const R3Point& position = pointset->PointPosition(i);
      const R3Vector& normal = pointset->PointNormal(i);
      RNScalar value = pointset->PointValue(i);
      if (show_colors) glColor3d(-10*value, 0, 10*value);
      else R3LoadNormal(normal);
      R3LoadPoint(position);
    }
    glEnd();
    glPointSize(1);
  }
  
  // Draw normals
  if (show_normals) {
    glDisable(GL_LIGHTING);
    glColor3d(0.0, 1.0, 0.0);
    glBegin(GL_LINES);
    RNScalar d = 0.0025 * world_bbox.DiagonalRadius();
    for (int i = 0; i < pointset->NPoints(); i++) {
      const R3Point& position = pointset->PointPosition(i);
      const R3Vector& normal = pointset->PointNormal(i);
      R3LoadPoint(position);
      R3LoadPoint(position + d * normal);
    }
    glEnd();
  }
}



static void
DrawGrid(R3Grid *grid, int model_index)
{
  // Push transformation
  grid->GridToWorldTransformation().Push();

  // Draw grid isosurface faces
  if (show_grid_faces) {
    glEnable(GL_LIGHTING);
    RNLoadRgb(0.8, 0.5, 0.2);
    grid->DrawIsoSurface(grid_thresholds[model_index]);
  }

  // Draw grid isosurface edges
  if (show_grid_edges) {
    glDisable(GL_LIGHTING);
    RNLoadRgb(0.0, 0.0, 0.0);
    glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
    grid->DrawIsoSurface(grid_thresholds[model_index]);
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
  }

  // Draw grid slices
  if (show_grid_slices[RN_X] || show_grid_slices[RN_Y] || show_grid_slices[RN_Z]) {
    glDisable(GL_LIGHTING);
    RNLoadRgb(1.0, 1.0, 1.0);
    if (show_grid_slices[RN_X]) grid->DrawSlice(RN_X, grid_slice_coords[RN_X]);
    if (show_grid_slices[RN_Y]) grid->DrawSlice(RN_Y, grid_slice_coords[RN_Y]);
    if (show_grid_slices[RN_Z]) grid->DrawSlice(RN_Z, grid_slice_coords[RN_Z]);
  }

  // Draw grid text
  if (show_text_labels) {
    char buffer[64];
    glDisable(GL_LIGHTING);
    RNLoadRgb(0.0, 0.0, 1.0);
    for (int i = 0; i < grid->XResolution(); i++) {
      for (int j = 0; j < grid->YResolution(); j++) {
        for (int k = 0; k < grid->ZResolution(); k++) {
          RNScalar value = grid->GridValue(i, j, k);
          if (value >= grid_thresholds[model_index]) continue;
          sprintf(buffer, "%.2g", value);
          DrawText(R3Point((RNScalar) i, (RNScalar) j, (RNScalar) k), buffer);
        }
      }
    }
  }

  // Draw grid threshold
  if (show_grid_threshold) {
    char buffer[128];
    sprintf(buffer, "%f", grid_thresholds[selected_model_index]);
    DrawText(R2Point(10,10), buffer);
  }

  // Pop transformation
  grid->GridToWorldTransformation().Pop();
}



static void
DrawBBox(void)
{
  // Draw bounding box
  glDisable(GL_LIGHTING);
  glColor3f(0.5, 0.5, 0.5);
  world_bbox.Outline();
}



static void
DrawAxes(void)
{
  // Draw axes
  RNScalar d = world_bbox.DiagonalRadius();
  glDisable(GL_LIGHTING);
  glLineWidth(3);
  R3BeginLine();
  glColor3f(1, 0, 0);
  R3LoadPoint(R3zero_point + 0.5*d * R3negx_vector);
  R3LoadPoint(R3zero_point + d * R3posx_vector);
  R3EndLine();
  R3BeginLine();
  glColor3f(0, 1, 0);
  R3LoadPoint(R3zero_point + 0.5*d * R3negy_vector);
  R3LoadPoint(R3zero_point + d * R3posy_vector);
  R3EndLine();
  R3BeginLine();
  glColor3f(0, 0, 1);
  R3LoadPoint(R3zero_point + 0.5*d * R3negz_vector);
  R3LoadPoint(R3zero_point + d * R3posz_vector);
  R3EndLine();
  glLineWidth(1);
}



////////////////////////////////////////////////////////////////////////
// GLUT interface functions
////////////////////////////////////////////////////////////////////////

static void
GLUTStop(void)
{
  // Destroy window 
  glutDestroyWindow(GLUTwindow);

  // Exit
  exit(0);
}



static void
GLUTRedraw(void)
{
  // Set viewing transformation
  viewer->Camera().Load();

  // Clear window 
  glClearColor(background_color[0], background_color[1], background_color[2], 1.0);
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

  // Set lights
  static GLfloat light0_position[] = { 3.0, 4.0, 5.0, 0.0 };
  glLightfv(GL_LIGHT0, GL_POSITION, light0_position);
  static GLfloat light1_position[] = { -3.0, -2.0, -3.0, 0.0 };
  glLightfv(GL_LIGHT1, GL_POSITION, light1_position);

  // Draw pointsets
  for (int m = 0; m < pointsets.NEntries(); m++) {
    if ((selected_model_index >= 0) && (m != selected_model_index)) continue;
    R3PointSet *pointset = pointsets[m];
    DrawPointSet(pointset);
  }
  
  // Draw scenes
  for (int m = 0; m < scenes.NEntries(); m++) {
    if ((selected_model_index >= 0) && (m != selected_model_index)) continue;
    R3Scene *scene = scenes[m];
    DrawScene(scene);
  }
  
  // Draw meshes
  for (int m = 0; m < meshes.NEntries(); m++) {
    if ((selected_model_index >= 0) && (m != selected_model_index)) continue;
    R3Mesh *mesh = meshes[m];
    DrawMesh(mesh);
  }
  
  // Draw grids
  for (int m = 0; m < grids.NEntries(); m++) {
    if ((selected_model_index >= 0) && (m != selected_model_index)) continue;
    R3Grid *grid = grids[m];
    DrawGrid(grid, m);
  }
  
  // Draw bbox
  if (show_bbox) {
    DrawBBox();
  }

  // Draw axes
  if (show_axes) {
    DrawAxes();
  }

  // Capture image and exit
  if (image_name) {
    R2Image image(GLUTwindow_width, GLUTwindow_height, 3);
    image.Capture();
    image.Write(image_name);
    GLUTStop();
  }

  // Swap buffers 
  glutSwapBuffers();
}    



static void
GLUTResize(int w, int h)
{
  // Resize window
  glViewport(0, 0, w, h);

  // Resize viewer viewport
  viewer->ResizeViewport(0, 0, w, h);

  // Remember window size 
  GLUTwindow_width = w;
  GLUTwindow_height = h;

  // Redraw
  glutPostRedisplay();
}



static void
GLUTMotion(int x, int y)
{
  // Invert y coordinate
  y = GLUTwindow_height - y;

  // Compute mouse movement
  int dx = x - GLUTmouse[0];
  int dy = y - GLUTmouse[1];
  
  // World in hand navigation 
  if (GLUTbutton[0]) viewer->RotateWorld(1.0, world_origin, x, y, dx, dy);
  else if (GLUTbutton[1]) viewer->ScaleWorld(1.0, world_origin, x, y, dx, dy);
  else if (GLUTbutton[2]) viewer->TranslateWorld(1.0, world_origin, x, y, dx, dy);
  if (GLUTbutton[0] || GLUTbutton[1] || GLUTbutton[2]) glutPostRedisplay();

  // Remember mouse position 
  GLUTmouse[0] = x;
  GLUTmouse[1] = y;
}



static void
GLUTMouse(int button, int state, int x, int y)
{
  // Invert y coordinate
  y = GLUTwindow_height - y;
  
  // Process mouse button event
  if ((button == GLUT_LEFT) && (state == GLUT_UP)) {
    // Check for double click  
    static RNBoolean double_click = FALSE;
    static RNTime last_mouse_down_time;
    double_click = (!double_click) && (last_mouse_down_time.Elapsed() < 0.4);
    last_mouse_down_time.Read();

    // Set world origin
    if (double_click) {
      R3Ray ray = viewer->WorldRay(x, y);
      RNScalar best_t = FLT_MAX;
      for (int m = 0; m < meshes.NEntries(); m++) {
        if ((selected_model_index >= 0) && (m != selected_model_index)) continue;
        R3Mesh *mesh = meshes[m];
        R3MeshIntersection intersection;
        if (mesh->Intersection(ray, &intersection)) {
          if (intersection.t < best_t) {
            world_origin = intersection.point;
            best_t = intersection.t;
          }
        }
      }
    }
  }

  // Remember button state 
  int b = (button == GLUT_LEFT_BUTTON) ? 0 : ((button == GLUT_MIDDLE_BUTTON) ? 1 : 2);
  GLUTbutton[b] = (state == GLUT_DOWN) ? 1 : 0;

  // Remember modifiers 
  GLUTmodifiers = glutGetModifiers();

  // Remember mouse position 
  GLUTmouse[0] = x;
  GLUTmouse[1] = y;

  // Redraw
  glutPostRedisplay();
}



static void
GLUTSpecial(int key, int x, int y)
{
  // Invert y coordinate
  y = GLUTwindow_height - y;

  // Process keyboard button event 
  switch (key) {
  case GLUT_KEY_DOWN:
    grid_thresholds[selected_model_index] -= 0.01 * grid_intervals[selected_model_index].Diameter();
    if (grid_thresholds[selected_model_index] < grid_intervals[selected_model_index].Min()) 
      grid_thresholds[selected_model_index] = grid_intervals[selected_model_index].Min();
    break;

  case GLUT_KEY_UP:
    grid_thresholds[selected_model_index] += 0.01 * grid_intervals[selected_model_index].Diameter();
    if (grid_thresholds[selected_model_index] > grid_intervals[selected_model_index].Max())
      grid_thresholds[selected_model_index] = grid_intervals[selected_model_index].Max();
    break;
      
  case GLUT_KEY_PAGE_UP:
    if (++selected_model_index > nmodels-1) selected_model_index = nmodels - 1;
    break;

  case GLUT_KEY_PAGE_DOWN:
    if (--selected_model_index < 0) selected_model_index = 0;
    break;
  }

  // Remember mouse position 
  GLUTmouse[0] = x;
  GLUTmouse[1] = y;

  // Remember modifiers 
  GLUTmodifiers = glutGetModifiers();

  // Redraw
  glutPostRedisplay();
}



static void
GLUTKeyboard(unsigned char key, int x, int y)
{
  // Process keyboard button event 
  switch (key) {
  case 'A':
  case 'a':
    show_axes = !show_axes;
    break;

  case 'B':
  case 'b':
    show_bbox = !show_bbox;
    break;

  case 'E':
  case 'e':
    show_scene_edges = !show_scene_edges;
    show_mesh_edges = !show_mesh_edges;
    break;

  case 'F':
  case 'f':
    show_scene_faces = !show_scene_faces;
    show_mesh_faces = !show_mesh_faces;
    break;

  case 'G':
  case 'g':
    show_grid_slices[0] = !show_grid_slices[0];
    show_grid_slices[1] = show_grid_slices[0];
    show_grid_slices[2] = show_grid_slices[0];
    break;

  case 'I':
    show_grid_faces = !show_grid_faces;
    break;

  case 'i':
    show_grid_edges = !show_grid_edges;
    break;

  case 'L':
  case 'l':
    show_categories = !show_categories;
    break;

  case 'N':
  case 'n':
    show_normals = !show_normals;
    break;

  case 'P':
  case 'p':
    show_points = !show_points;
    break;

  case 'R':
  case 'r':
    show_colors = !show_colors;
    break;

  case 'S':
  case 's':
    show_segments = !show_segments;
    break;

  case 'T':
    show_grid_threshold = !show_grid_threshold;
    break;

  case 't':
    show_text_labels = !show_text_labels;
    break;

  case 'V':
  case 'v':
    show_mesh_vertices = !show_mesh_vertices;
    break;

  case 'X':
    if (selected_model_index < grids.NEntries()) {
      grid_slice_coords[RN_X]++;
      if (grid_slice_coords[RN_X] >= grids[selected_model_index]->XResolution()) 
        grid_slice_coords[RN_X] = grids[selected_model_index]->XResolution() - 1;
    }
    break;

  case 'x':
    if (selected_model_index < grids.NEntries()) {
      grid_slice_coords[RN_X]--;
      if (grid_slice_coords[RN_X] < 0)
        grid_slice_coords[RN_X] = 0;
    }
    break;

  case 'Y':
    if (selected_model_index < grids.NEntries()) {
      grid_slice_coords[RN_Y]++;
      if (grid_slice_coords[RN_Y] >= grids[selected_model_index]->YResolution()) 
        grid_slice_coords[RN_Y] = grids[selected_model_index]->YResolution() - 1;
    }
    break;

  case 'y':
    if (selected_model_index < grids.NEntries()) {
      grid_slice_coords[RN_Y]--;
      if (grid_slice_coords[RN_Y] < 0)
        grid_slice_coords[RN_Y] = 0;
    }
    break;

  case 'Z':
    if (selected_model_index < grids.NEntries()) {
      grid_slice_coords[RN_Z]++;
      if (grid_slice_coords[RN_Z] >= grids[selected_model_index]->ZResolution()) 
        grid_slice_coords[RN_Z] = grids[selected_model_index]->ZResolution() - 1;
    }
    break;

  case 'z':
    if (selected_model_index < grids.NEntries()) {
      grid_slice_coords[RN_Z]--;
      if (grid_slice_coords[RN_Z] < 0)
        grid_slice_coords[RN_Z] = 0;
    }
    break;

  case ' ': {
    // Print camera
    const R3Camera& camera = viewer->Camera();
    printf("#camera  %g %g %g  %g %g %g  %g %g %g  %g \n",
           camera.Origin().X(), camera.Origin().Y(), camera.Origin().Z(),
           camera.Towards().X(), camera.Towards().Y(), camera.Towards().Z(),
           camera.Up().X(), camera.Up().Y(), camera.Up().Z(),
           camera.YFOV());
    break; }
      
  case 27: // ESCAPE
    GLUTStop();
    break;
  }

  // Remember mouse position 
  GLUTmouse[0] = x;
  GLUTmouse[1] = GLUTwindow_height - y;

  // Remember modifiers 
  GLUTmodifiers = glutGetModifiers();

  // Redraw
  glutPostRedisplay();  
}




static void
GLUTInit(int *argc, char **argv)
{
  // Open window 
  glutInit(argc, argv);
  glutInitWindowPosition(100, 100);
  glutInitWindowSize(GLUTwindow_width, GLUTwindow_height);
  glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH); // | GLUT_STENCIL
  GLUTwindow = glutCreateWindow("OpenGL Viewer");

  // Initialize background color 
  glClearColor(200.0/255.0, 200.0/255.0, 200.0/255.0, 1.0);

  // Initialize lights
  static GLfloat lmodel_ambient[] = { 0.2, 0.2, 0.2, 1.0 };
  glLightModelfv(GL_LIGHT_MODEL_AMBIENT, lmodel_ambient);
  glLightModeli(GL_LIGHT_MODEL_LOCAL_VIEWER, GL_TRUE);
  static GLfloat light0_diffuse[] = { 1.0, 1.0, 1.0, 1.0 };
  glLightfv(GL_LIGHT0, GL_DIFFUSE, light0_diffuse);
  glEnable(GL_LIGHT0);
  static GLfloat light1_diffuse[] = { 0.5, 0.5, 0.5, 1.0 };
  glLightfv(GL_LIGHT1, GL_DIFFUSE, light1_diffuse);
  glEnable(GL_LIGHT1);
  glEnable(GL_NORMALIZE);
  glEnable(GL_LIGHTING);

  // Initialize graphics modes 
  glEnable(GL_DEPTH_TEST);

  // Initialize GLUT callback functions 
  glutDisplayFunc(GLUTRedraw);
  glutReshapeFunc(GLUTResize);
  glutKeyboardFunc(GLUTKeyboard);
  glutSpecialFunc(GLUTSpecial);
  glutMouseFunc(GLUTMouse);
  glutMotionFunc(GLUTMotion);

  // Initialize font
#if (RN_OS == RN_WINDOWSNT)
  int font = glGenLists(256);
  wglUseFontBitmaps(wglGetCurrentDC(), 0, 256, font); 
  glListBase(font);
#endif
}



static R3Viewer *
CreateViewer(void)
{
  // Start statistics
  RNTime start_time;
  start_time.Read();

  // Get pdb bounding box
  assert(!world_bbox.IsEmpty());
  RNLength r = world_bbox.DiagonalRadius();
  assert((r > 0.0) && RNIsFinite(r));

  // Setup camera view looking down the Z axis
  if (!initial_camera) initial_camera_origin = world_bbox.Centroid() - initial_camera_towards * (2.5 * r);;
  R3Camera camera(initial_camera_origin, initial_camera_towards, initial_camera_up, 0.4, 0.4, 0.1 * r, 1000.0 * r);
  R2Viewport viewport(0, 0, GLUTwindow_width, GLUTwindow_height);
  R3Viewer *viewer = new R3Viewer(camera, viewport);

  // Print statistics
  if (print_verbose) {
    printf("Created viewer ...\n");
    printf("  Time = %.2f seconds\n", start_time.Elapsed());
    printf("  Origin = %g %g %g\n", camera.Origin().X(), camera.Origin().Y(), camera.Origin().Z());
    printf("  Towards = %g %g %g\n", camera.Towards().X(), camera.Towards().Y(), camera.Towards().Z());
    printf("  Up = %g %g %g\n", camera.Up().X(), camera.Up().Y(), camera.Up().Z());
    printf("  Fov = %g %g\n", camera.XFOV(), camera.YFOV());
    printf("  Near = %g\n", camera.Near());
    printf("  Far = %g\n", camera.Far());
    fflush(stdout);
  }

  // Return viewer
  return viewer;
}



void GLUTMainLoop(void)
{
 // Set world bbbox
  world_bbox = R3null_box;
  for (int i = 0; i < scenes.NEntries(); i++) world_bbox.Union(scenes[i]->BBox());
  for (int i = 0; i < meshes.NEntries(); i++) world_bbox.Union(meshes[i]->BBox());
  for (int i = 0; i < pointsets.NEntries(); i++) world_bbox.Union(pointsets[i]->BBox());
  for (int i = 0; i < grids.NEntries(); i++) world_bbox.Union(grids[i]->WorldBox());

  // Set world origin
  world_origin = world_bbox.Centroid();

  // Set slice coords to middle of grid
  if (grids.NEntries() > 0) {
    grid_slice_coords[RN_X] = grids[0]->XResolution()/2;
    grid_slice_coords[RN_Y] = grids[0]->YResolution()/2;
    grid_slice_coords[RN_Z] = grids[0]->ZResolution()/2;
  }

  // Create viewer
  viewer = CreateViewer();
  if (!viewer) exit(-1);

   // Run main loop -- never returns 
  glutMainLoop();
}



static int 
ParseArgs(int argc, char **argv)
{
  // Parse arguments
  argc--; argv++;
  while (argc > 0) {
    if ((*argv)[0] == '-') {
      if (!strcmp(*argv, "-v")) print_verbose = 1; 
      else if (!strcmp(*argv, "-grid_threshold")) {
        argc--; argv++; 
        if (grids.NEntries() > 0) {
          grid_thresholds[grids.NEntries()-1] = atof(*argv);
        }
      }
      else if (!strcmp(*argv, "-camera")) { 
        RNCoord x, y, z, tx, ty, tz, ux, uy, uz;
        argv++; argc--; x = atof(*argv); 
        argv++; argc--; y = atof(*argv); 
        argv++; argc--; z = atof(*argv); 
        argv++; argc--; tx = atof(*argv); 
        argv++; argc--; ty = atof(*argv); 
        argv++; argc--; tz = atof(*argv); 
        argv++; argc--; ux = atof(*argv); 
        argv++; argc--; uy = atof(*argv); 
        argv++; argc--; uz = atof(*argv); 
        initial_camera_origin = R3Point(x, y, z);
        initial_camera_towards.Reset(tx, ty, tz);
        initial_camera_up.Reset(ux, uy, uz);
        initial_camera = TRUE;
      }
      else if (!strcmp(*argv, "-background")) { 
        argc--; argv++; background_color[0] = atof(*argv); 
        argc--; argv++; background_color[1] = atof(*argv); 
        argc--; argv++; background_color[2] = atof(*argv); 
      }
      else if (!strcmp(*argv, "-image")) { 
        argc--; argv++; image_name = *argv; 
      }
      else { 
        RNFail("Invalid program argument: %s", *argv); 
        exit(1); 
      }
      argv++; argc--;
    }
    else {
      char *ext = strrchr(*argv, '.');
      if (ext && (!strcmp(ext, ".obj") || !strcmp(ext, ".scn") || !strcmp(ext, ".json"))) {
        R3Scene *scene = ReadScene(*argv);
        if (!scene) return 0;
        scenes.Insert(scene);
        if (scenes.NEntries() > nmodels){
          nmodels = scenes.NEntries();
        }
      }
      else if (ext && (!strcmp(ext, ".off") || !strcmp(ext, ".ply") || !strcmp(ext, ".wrl"))) {
        R3Mesh *mesh = ReadMesh(*argv);
        if (!mesh) return 0;
        meshes.Insert(mesh);
        if (meshes.NEntries() > nmodels){
          nmodels = meshes.NEntries();
        }
      }
      else if (ext && (!strcmp(ext, ".xyzn") || !strcmp(ext, ".pts") || !strcmp(ext, ".sdf"))) {
        R3PointSet *pointset = ReadPointSet(*argv);
        if (!pointset) return 0;
        pointsets.Insert(pointset);
        if (pointsets.NEntries() > nmodels){
          nmodels = pointsets.NEntries();
        }
      }
      else if (ext && (!strcmp(ext, ".grd"))) {
        R3Grid *grid = ReadGrid(*argv);
        if (!grid) return 0;
        grid_intervals[grids.NEntries()] = grid->Range();
        grid_thresholds[grids.NEntries()] = grid_intervals[grids.NEntries()].Mid();
        grids.Insert(grid);
        if (grids.NEntries() > nmodels){
          nmodels = grids.NEntries();
        }
      }
      else { 
        RNFail("Invalid program argument: %s", *argv); 
        exit(1); 
      }
      argv++; argc--;
    }
  }

  // Check inputs
  if (nmodels == 0) {
    RNFail("Usage: view <scenes> <meshes> <grids> [options]\n");
    return 0;
  }

  // Return OK status 
  return 1;
}



int 
main(int argc, char **argv)
{
  // Initialize GLUT
  GLUTInit(&argc, argv);

  // Parse program arguments
  if (!ParseArgs(argc, argv)) exit(-1);

  // Run GLUT interface
  GLUTMainLoop();

  // Return success 
  return 0;
}



