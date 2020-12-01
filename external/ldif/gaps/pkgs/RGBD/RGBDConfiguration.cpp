////////////////////////////////////////////////////////////////////////
// Source file for RGBDConfiguration class
////////////////////////////////////////////////////////////////////////



////////////////////////////////////////////////////////////////////////
// Include files
////////////////////////////////////////////////////////////////////////

#include "RGBD.h"




////////////////////////////////////////////////////////////////////////
// Namespace
////////////////////////////////////////////////////////////////////////

namespace gaps {



////////////////////////////////////////////////////////////////////////
// Constructors/destructors
////////////////////////////////////////////////////////////////////////

RGBDConfiguration::
RGBDConfiguration(void)
  : images(),
    surfaces(),
    name(NULL),
    color_directory(NULL),
    depth_directory(NULL),
    category_directory(NULL),
    instance_directory(NULL),
    texture_directory(NULL),
    dataset_format(NULL),
    world_bbox(FLT_MAX, FLT_MAX, FLT_MAX, -FLT_MAX, -FLT_MAX, -FLT_MAX)
{
}



RGBDConfiguration::
~RGBDConfiguration(void)
{
  // Delete everything
  while (NImages() > 0) delete Image(NImages()-1);
  while (NSurfaces() > 0) delete Surface(NSurfaces()-1);

  // Delete name
  if (name) free(name);

  // Delete directory names
  if (color_directory) free(color_directory);
  if (depth_directory) free(depth_directory);
  if (category_directory) free(category_directory);
  if (instance_directory) free(instance_directory);
  if (texture_directory) free(texture_directory);

  // Delete dataset format
  if (dataset_format) free(dataset_format);
}



////////////////////////////////////////////////////////////////////////
// Property functions
////////////////////////////////////////////////////////////////////////

R3Point RGBDConfiguration::
WorldCentroid(void) const
{
  // Compute weighed sum of surface centroids
  RNArea weight = 0;
  R3Point position = R3zero_point;
  for (int i = 0; i < NSurfaces(); i++) {
    RGBDSurface *surface = Surface(i);
    position += surface->WorldArea() * surface->WorldCentroid();
    weight += surface->WorldArea();
  }

  // Return weighted average
  if (RNIsZero(weight)) return R3zero_point;
  return position / weight;
}



////////////////////////////////////////////////////////////////////////
// Image and surface manipulation functions
////////////////////////////////////////////////////////////////////////

void RGBDConfiguration::
InsertImage(RGBDImage *image)
{
  // Just checking
  assert(image);
  assert(image->configuration == NULL);
  assert(image->configuration_index == -1);

  // Insert image into configuration
  image->configuration = this;
  image->configuration_index = images.NEntries();
  images.Insert(image);

  // Update bounding box
  world_bbox.Union(image->WorldBBox());
}



void RGBDConfiguration::
RemoveImage(RGBDImage *image)
{
  // Just checking
  assert(image);
  assert(image->configuration == this);
  assert(image->configuration_index >= 0);
  assert(images.Kth(image->configuration_index) == image);

  // Remove image from configuration
  RNArrayEntry *entry = images.KthEntry(image->configuration_index);
  RGBDImage *tail = images.Tail();
  tail->configuration_index = image->configuration_index;
  images.EntryContents(entry) = tail;
  images.RemoveTail();
  image->configuration_index = -1;
  image->configuration = NULL;
}



void RGBDConfiguration::
InsertSurface(RGBDSurface *surface)
{
  // Just checking
  assert(surface);
  assert(surface->configuration == NULL);
  assert(surface->configuration_index == -1);

  // Insert surface into configuration
  surface->configuration = this;
  surface->configuration_index = surfaces.NEntries();
  surfaces.Insert(surface);

  // Update bounding box
  world_bbox.Union(surface->WorldBBox());
}



void RGBDConfiguration::
RemoveSurface(RGBDSurface *surface)
{
  // Just checking
  assert(surface);
  assert(surface->configuration == this);
  assert(surface->configuration_index >= 0);
  assert(surfaces.Kth(surface->configuration_index) == surface);

  // Remove surface from configuration
  RNArrayEntry *entry = surfaces.KthEntry(surface->configuration_index);
  RGBDSurface *tail = surfaces.Tail();
  tail->configuration_index = surface->configuration_index;
  surfaces.EntryContents(entry) = tail;
  surfaces.RemoveTail();
  surface->configuration_index = -1;
  surface->configuration = NULL;
}



void RGBDConfiguration::
Transform(const R3Transformation& transformation)
{
  // Transform images
  for (int i = 0; i < NImages(); i++) {
    RGBDImage *image = Image(i);
    image->Transform(transformation);
  }
  
  // Transform surfaces
  for (int i = 0; i < NSurfaces(); i++) {
    RGBDSurface *surface = Surface(i);
    surface->Transform(transformation);
  }
}



////////////////////////////////////////////////////////////////////////
// Directory name manipulation functions
////////////////////////////////////////////////////////////////////////

void RGBDConfiguration::
SetName(const char *name)
{
  // Set directory name
  if (this->name) free(this->name);
  if (name && strcmp(name, "-")) this->name = RNStrdup(name);
  else this->name = NULL;
}



void RGBDConfiguration::
SetColorDirectory(const char *directory)
{
  // Set directory name
  if (color_directory) free(color_directory);
  if (directory && strcmp(directory, "-")) color_directory = RNStrdup(directory);
  else color_directory = NULL;
}



void RGBDConfiguration::
SetDepthDirectory(const char *directory)
{
  // Set directory name
  if (depth_directory) free(depth_directory);
  if (directory && strcmp(directory, "-")) depth_directory = RNStrdup(directory);
  else depth_directory = NULL;
}



void RGBDConfiguration::
SetCategoryDirectory(const char *directory)
{
  // Set directory name
  if (category_directory) free(category_directory);
  if (directory && strcmp(directory, "-")) category_directory = RNStrdup(directory);
  else category_directory = NULL;
}



void RGBDConfiguration::
SetInstanceDirectory(const char *directory)
{
  // Set directory name
  if (instance_directory) free(instance_directory);
  if (directory && strcmp(directory, "-")) instance_directory = RNStrdup(directory);
  else instance_directory = NULL;
}



void RGBDConfiguration::
SetTextureDirectory(const char *directory)
{
  // Set directory name
  if (texture_directory) free(texture_directory);
  if (directory && strcmp(directory, "-")) texture_directory = RNStrdup(directory);
  else texture_directory = NULL;
}



void RGBDConfiguration::
SetDatasetFormat(const char *format)
{
  // Set format name
  if (dataset_format) free(dataset_format);
  if (format && strcmp(format, "-")) dataset_format = RNStrdup(format);
  else dataset_format = NULL;
}



////////////////////////////////////////////////////////////////////////
// File input/output functions
////////////////////////////////////////////////////////////////////////

int RGBDConfiguration::
ReadFile(const char *filename, int read_every_kth_image)
{
  // Parse input filename extension
  const char *extension;
  if (!(extension = strrchr(filename, '.'))) {
    printf("Filename %s has no extension (e.g., .conf)\n", filename);
    return 0;
  }

  // Read file of appropriate type
  if (!strncmp(extension, ".conf", 5)) {
    if (!ReadConfigurationFile(filename, read_every_kth_image)) return 0;
  }
  else {
    RNFail("Unable to read file %s (unrecognized extension: %s)\n", filename, extension);
    return 0;
  }

  // Return success
  return 1;
}


int RGBDConfiguration::
WriteFile(const char *filename, int write_every_kth_image) const
{
  // Parse input filename extension
  const char *extension;
  if (!(extension = strrchr(filename, '.'))) {
    printf("Filename %s has no extension (e.g., .conf)\n", filename);
    return 0;
  }

  // Write file of appropriate type
  if (!strncmp(extension, ".conf", 5)) {
    if (!WriteConfigurationFile(filename, write_every_kth_image)) return 0;
  }
  else if (!strncmp(extension, ".obj", 4)) {
    if (!WriteObjFile(filename)) return 0;
  }
  else {
    RNFail("Unable to write file %s (unrecognized extension: %s)\n", filename, extension);
    return 0;
  }

  // Return success
  return 1;
}


////////////////////////////////////////////////////////////////////////
// Configuration file input/output functions
////////////////////////////////////////////////////////////////////////

int RGBDConfiguration::
ReadConfigurationFile(const char *filename, int read_every_kth_image)
{
  // Open file
  FILE *fp = fopen(filename, "r");
  if (!fp) {
    RNFail("Unable to open configuration file %s\n", filename);
    return 0;
  }

  // Read file
  if (!ReadConfigurationStream(fp, read_every_kth_image)) return 0;
  
  // Close file
  fclose(fp);

  // Return success
  return 1;
}



int RGBDConfiguration::
WriteConfigurationFile(const char *filename, int write_every_kth_image) const
{
#if 0
  // Create directories
  char cmd[4096];
  if (color_directory) { sprintf(cmd, "mkdir -p %s", color_directory); system(cmd); }
  if (depth_directory) { sprintf(cmd, "mkdir -p %s", depth_directory); system(cmd); }
  if (category_directory) { sprintf(cmd, "mkdir -p %s", category_directory); system(cmd); }
  if (instance_directory) { sprintf(cmd, "mkdir -p %s", instance_directory); system(cmd); }
  if (texture_directory) { sprintf(cmd, "mkdir -p %s", texture_directory); system(cmd); }
#endif

  // Open file
  FILE *fp = fopen(filename, "w");
  if (!fp) {
    RNFail("Unable to open configuration file %s\n", filename);
    return 0;
  }

  // Write file
  if (!WriteConfigurationStream(fp, write_every_kth_image)) return 0;

  // Close file
  fclose(fp);

  // Return success
  return 1;
}



int RGBDConfiguration::
ReadConfigurationStream(FILE *fp, int read_every_kth_image)
{
  // Initialize stuff
  int image_count = 0;
  R3Matrix intrinsics_matrix = R3identity_matrix;
  int image_width = 0;
  int image_height = 0;
  
  // Parse file
  char buffer[4096];
  int line_number = 0;
  while (fgets(buffer, 4096, fp)) {
    char cmd[1024];
    line_number++;
    if (sscanf(buffer, "%s", cmd) != (unsigned int) 1) continue;
    if (cmd[0] == '#') continue;

    // Check cmd
    if (!strcmp(cmd, "dataset")) {
      // Parse dataset format
      char name[1024] = { '\0' };
      if (sscanf(buffer, "%s%s", cmd, name) != (unsigned int) 2) {
        RNFail("Error parsing line %d of configuration file\n", line_number);
        return 0;
      }

      // Set dataset format
      SetDatasetFormat(name);
    }
    else if (!strcmp(cmd, "sequence") || !strcmp(cmd, "name")) {
      // Parse name
      char name[1024] = { '\0' };
      if (sscanf(buffer, "%s%s", cmd, name) != (unsigned int) 2) {
        RNFail("Error parsing line %d of configuration file\n", line_number);
        return 0;
      }

      // Set dataset format
      SetName(name);
    }
    else if (!strcmp(cmd, "depth_resolution")) {
      // Parse dimensions
      if (sscanf(buffer, "%s%d%d", cmd, &image_width, &image_height) != (unsigned int) 3) {
        RNFail("Error parsing line %d of configuration file\n", line_number);
        return 0;
      }
    }
    else if (!strcmp(cmd, "color_resolution")) {
      // Parse dimensions
      if ((image_width == 0) || (image_height == 0)) {
        if (sscanf(buffer, "%s%d%d", cmd, &image_width, &image_height) != (unsigned int) 3) {
          RNFail("Error parsing line %d of configuration file\n", line_number);
          return 0;
        }
      }
    }
    else if (!strcmp(cmd, "depth_directory") || !strcmp(cmd, "color_directory") ||
             !strcmp(cmd, "category_directory") || !strcmp(cmd, "instance_directory") ||
             !strcmp(cmd, "texture_directory") || !strcmp(cmd, "image_directory")) {
      // Parse directory name
      char dirname[1024];
      if (sscanf(buffer, "%s%s", cmd, dirname) != (unsigned int) 2) {
        RNFail("Error parsing line %d of configuration file\n", line_number);
        return 0;
      }

      // Assign directory name (this leaks memory, but who cares)
      if (!strcmp(cmd, "color_directory")) SetColorDirectory(dirname);
      else if (!strcmp(cmd, "depth_directory")) SetDepthDirectory(dirname);
      else if (!strcmp(cmd, "category_directory")) SetCategoryDirectory(dirname);
      else if (!strcmp(cmd, "instance_directory")) SetInstanceDirectory(dirname);
      else if (!strcmp(cmd, "image_directory")) SetColorDirectory(dirname);
      else if (!strcmp(cmd, "texture_directory")) SetTextureDirectory(dirname);
    }
    else if (!strcmp(cmd, "intrinsics") || !strcmp(cmd, "depth_intrinsics")) { // || !strcmp(cmd, "color_intrinsics")) {
      // Parse intrinsics filename
      char intrinsics_filename[2048];
      if (sscanf(buffer, "%s%s", cmd, intrinsics_filename) != (unsigned int) 2) {
        RNFail("Error parsing line %d of configuration file\n", line_number);
        return 0;
      }

      // Open intrinsics file
      FILE *intrinsics_fp = fopen(intrinsics_filename, "r");
      if (!intrinsics_fp) {
        RNFail("Unable to open intrinsics file %s\n", intrinsics_filename);
        return 0;
      }

      // Read intrinsics file
      if (dataset_format && !strcmp(dataset_format, "matterport")) {
        // Read matterport intrinsics
        double width, height, fx, fy, cx, cy, k1, k2, p1, p2, k3;
        if (fscanf(intrinsics_fp, "%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf", &width, &height,
          &fx, &fy, &cx, &cy, &k1, &k2, &p1, &p2, &k3) != (unsigned int) 11) {
          RNFail("Unable to read Matterport intrinsics matrix.\n");
          return 0;
        }

        // Assign intrinsics matrix
        intrinsics_matrix = R3Matrix(fx, 0, cx,   0, fy, cy,   0, 0, 1);
      }
      else {
        // Read matrix
        RNScalar m[9];
        if (fscanf(intrinsics_fp, "%lf%lf%lf%lf%lf%lf%lf%lf%lf", &m[0], &m[1], &m[2], &m[3], &m[4], &m[5], &m[6], &m[7], &m[8]) != (unsigned int) 9) {
          RNFail("Unable to read intrinsics file %s\n", intrinsics_filename);
          return 0;
        }

        // Assign intrinsics matrix
        intrinsics_matrix = R3Matrix(m);
      }

      // Close intrinsics file
      fclose(intrinsics_fp);
    }
    else if (!strcmp(cmd, "intrinsics_matrix") || !strcmp(cmd, "depth_intrinsics_matrix")) { // || !strcmp(cmd, "color_intrinsics_matrix")) {
      // Parse matrix
      double m[9];
      if (sscanf(buffer, "%s%lf%lf%lf%lf%lf%lf%lf%lf%lf", cmd, &m[0], &m[1], &m[2], &m[3], &m[4], &m[5], &m[6], &m[7], &m[8]) != (unsigned int) 10) {
        RNFail("Error parsing line %d of configuration file\n", line_number);
        return 0;
      }

      // Assign intrinsics matrix
      intrinsics_matrix = R3Matrix(m);
    }
    else if (!strcmp(cmd, "scan") || !strcmp(cmd, "image")) {
      // Update/check image count
      if ((read_every_kth_image > 1) && ((++image_count % read_every_kth_image) != 1)) continue;

      // Check intrinsics matrix
      if (intrinsics_matrix.IsIdentity()) {
        RNFail("Unable to process scan without prior setting of intrinsics matrix in configuration file\n");
        return 0;
      }
      
      // Parse image names and alignment transformation
      RNScalar m[16];
      char depth_filename[2048], color_filename[2048];
      if (sscanf(buffer, "%s%s%s%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf", cmd, 
         depth_filename, color_filename,
         &m[0], &m[1], &m[2], &m[3], &m[4], &m[5], &m[6], &m[7], 
         &m[8], &m[9], &m[10], &m[11], &m[12], &m[13], &m[14], &m[15]) != (unsigned int) 19) {
        RNFail("Error parsing line %d of configuration file\n", line_number);
        return 0;
      }

      // Create RGBD image
      // RGBDImage *image = new RGBDImage(color_filename, depth_filename, intrinsics_matrix, R4Matrix(m), image_width, image_height);
      RGBDImage *image = AllocateImage();
      image->SetNPixels(image_width, image_height);
      image->SetColorFilename(color_filename);
      image->SetDepthFilename(depth_filename);
      image->SetIntrinsics(intrinsics_matrix);
      image->SetCameraToWorld(R3Affine(R4Matrix(m), 0));
      InsertImage(image);
    }
    else if (!strcmp(cmd, "labeled_image")) {
      // Update/check image count
      if ((read_every_kth_image > 1) && ((++image_count % read_every_kth_image) != 1)) continue;

      // Check intrinsics matrix
      if (intrinsics_matrix.IsIdentity()) {
        RNFail("Unable to process scan without prior setting of intrinsics matrix in configuration file\n");
        return 0;
      }
      
      // Parse image names and alignment transformation
      RNScalar m[16];
      char depth_filename[2048], color_filename[2048], category_filename[2048], instance_filename[2048];
      if (sscanf(buffer, "%s%s%s%s%s%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf", cmd, 
         depth_filename, color_filename, category_filename, instance_filename,
         &m[0], &m[1], &m[2], &m[3], &m[4], &m[5], &m[6], &m[7], 
         &m[8], &m[9], &m[10], &m[11], &m[12], &m[13], &m[14], &m[15]) != (unsigned int) 21) {
        RNFail("Error parsing line %d of configuration file\n", line_number);
        return 0;
      }

      // Create RGBD image
      // RGBDImage *image = new RGBDImage(color_filename, depth_filename, intrinsics_matrix, R4Matrix(m), image_width, image_height);
      RGBDImage *image = AllocateImage();
      image->SetNPixels(image_width, image_height);
      image->SetColorFilename(color_filename);
      image->SetDepthFilename(depth_filename);
      image->SetCategoryFilename(category_filename);
      image->SetInstanceFilename(instance_filename);
      image->SetIntrinsics(intrinsics_matrix);
      image->SetCameraToWorld(R3Affine(R4Matrix(m), 0));
      InsertImage(image);
    }
    else if (!strcmp(cmd, "frame")) {
      // Update/check image count
      if ((read_every_kth_image > 1) && ((++image_count % read_every_kth_image) != 1)) continue;

      // Check intrinsics matrix
      if (intrinsics_matrix.IsIdentity()) {
        RNFail("Unable to process scan without prior setting of intrinsics matrix in configuration file\n");
        return 0;
      }
      
      // Parse image names and alignment transformation
      RNScalar m[16];
      double depth_timestamp, color_timestamp;
      char depth_filename[2048], color_filename[2048];
      if (sscanf(buffer, "%s%s%lf%s%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf", cmd, 
         depth_filename, &depth_timestamp, color_filename, &color_timestamp,
         &m[0], &m[1], &m[2], &m[3], &m[4], &m[5], &m[6], &m[7], 
         &m[8], &m[9], &m[10], &m[11], &m[12], &m[13], &m[14], &m[15]) != (unsigned int) 21) {
        RNFail("Error parsing line %d of configuration file\n", line_number);
        return 0;
      }

      // Create RGBD image
      // RGBDImage *image = new RGBDImage(color_filename, depth_filename, intrinsics_matrix, R4Matrix(m), image_width, image_height);
      RGBDImage *image = AllocateImage();
      image->SetNPixels(image_width, image_height);
      image->SetColorFilename(color_filename);
      image->SetDepthFilename(depth_filename);
      image->SetIntrinsics(intrinsics_matrix);
      image->SetCameraToWorld(R3Affine(R4Matrix(m), 0));
      InsertImage(image);
    }
    else if (!strcmp(cmd, "rectangle")) {
      // Parse surface
      char texture_filename[2048];
      RNScalar pixel_spacing, c[3], n[3], u[3], r[2];
      if (sscanf(buffer, "%s%s%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf", cmd, texture_filename, &pixel_spacing, 
         &c[0], &c[1], &c[2], &n[0], &n[1], &n[2], &u[0], &u[1], &u[2], &r[0], &r[1]) != (unsigned int) 14) {
        RNFail("Error parsing line %d of configuration file\n", line_number);
        return 0;
      }

      // Create RGBD surface
      R3Point center(c);
      R3Vector axis1(u); axis1.Normalize();
      R3Vector normal(n); normal.Normalize();
      R3Vector axis0 = axis1 % normal; axis0.Normalize();
      // R3Mesh *mesh = new R3Mesh();
      // R3MeshVertex *v00 = mesh->CreateVertex(center - r[0]*axis0 - r[1]*axis1, normal, RNblack_rgb, R2Point(0,0));
      // R3MeshVertex *v10 = mesh->CreateVertex(center + r[0]*axis0 - r[1]*axis1, normal, RNblack_rgb, R2Point(r[0],0));
      // R3MeshVertex *v01 = mesh->CreateVertex(center - r[0]*axis0 + r[1]*axis1, normal, RNblack_rgb, R2Point(0,r[1]));
      // R3MeshVertex *v11 = mesh->CreateVertex(center + r[0]*axis0 + r[1]*axis1, normal, RNblack_rgb, R2Point(r[0],r[1]));
      // mesh->CreateFace(v00, v10, v11);
      // mesh->CreateFace(v00, v11, v01);
      // RGBDSurface *surface = new RGBDSurface(texture_filename, mesh, pixel_spacing);
      R3Rectangle *rectangle = new R3Rectangle(center, axis0, axis1, r[0], r[1]);
      RGBDSurface *surface = new RGBDSurface(texture_filename, rectangle, pixel_spacing);
      InsertSurface(surface);
    }
    else if (!strcmp(cmd, "mesh")) {
      // Parse surface
      char texture_filename[2048], mesh_filename[2048];
      RNScalar pixel_spacing;
      if (sscanf(buffer, "%s%s%s%lf", cmd, texture_filename, mesh_filename, &pixel_spacing) != (unsigned int) 4) {
        RNFail("Error parsing line %d of configuration file\n", line_number);
        return 0;
      }

      // Read mesh
      R3Mesh *mesh = new R3Mesh();
      if (!mesh->ReadFile(mesh_filename)) {
        RNFail("Unable to read mesh file %s at line %d of configuration file\n", mesh_filename, line_number);
        return 0;
      }

      // Create surface      
      RGBDSurface *surface = new RGBDSurface(texture_filename, mesh, pixel_spacing);
      surface->SetMeshFilename(mesh_filename);
      InsertSurface(surface);
    }
  }

  // Return success
  return 1;
}



int RGBDConfiguration::
WriteConfigurationStream(FILE *fp, int write_every_kth_image) const
{
  // Write header
  fprintf(fp, "n_images %d\n", NImages());
  if (Name()) fprintf(fp, "sequence %s\n", Name());
  if (DatasetFormat()) fprintf(fp, "dataset %s\n", DatasetFormat());
  if (color_directory) fprintf(fp, "color_directory %s\n", color_directory);
  if (depth_directory) fprintf(fp, "depth_directory %s\n", depth_directory);
  if (category_directory) fprintf(fp, "category_directory %s\n", category_directory);
  if (instance_directory) fprintf(fp, "instance_directory %s\n", instance_directory);
  if (texture_directory) fprintf(fp, "texture_directory %s\n", texture_directory);
  
  // Write blank line
  fprintf(fp, "\n");

  // Write surfaces
  for (int i = 0; i < NSurfaces(); i++) {
    RGBDSurface *surface = Surface(i);
    // surface->WriteChannels();
    if (surface->rectangle) {
      R3Point centroid = surface->rectangle->Centroid();
      R3Vector normal = surface->rectangle->Normal();
      R3Vector up = surface->rectangle->Axis(1);
      RNLength rx = surface->rectangle->Radius(0);
      RNLength ry = surface->rectangle->Radius(1);
      const char *texture_filename = (surface->TextureFilename()) ? surface->TextureFilename() : "-";
      fprintf(fp, "rectangle  %s  %g   %g %g %g   %g %g %g   %g %g %g   %g %g\n",
        texture_filename, surface->WorldTexelSpacing(),
        centroid.X(), centroid.Y(), centroid.Z(),
        normal.X(), normal.Y(), normal.Z(),
        up.X(), up.Y(), up.Z(), rx, ry);
    }
    else if (surface->mesh) {
      const char *texture_filename = (surface->TextureFilename()) ? surface->TextureFilename() : "-";
      const char *mesh_filename = (surface->MeshFilename()) ? surface->MeshFilename() : "-";
      fprintf(fp, "mesh  %s  %s  %g\n", texture_filename, mesh_filename, surface->WorldTexelSpacing());
    }
  }
  
  // Write blank line
  fprintf(fp, "\n");

  // Write images
  int image_width = 0;
  int image_height = 0;
  R3Matrix intrinsics_matrix = R3identity_matrix;
  for (int i = 0; i < NImages(); i++) {
    RGBDImage *image = Image(i);
    if ((write_every_kth_image > 1) && ((i % write_every_kth_image) != 0)) continue;

    // Write image dimensions
    if ((image->NPixels(RN_X) > 0) && (image->NPixels(RN_Y) > 0)) {
      if ((image->NPixels(RN_X) != image_width) || (image->NPixels(RN_Y) != image_height)) {
        fprintf(fp, "depth_resolution %d %d\n", image->NPixels(RN_X), image->NPixels(RN_Y));
        image_width = image->NPixels(RN_X);
        image_height = image->NPixels(RN_Y);
      }
    }

    // Write intrinsics matrix
    if (!image->Intrinsics().IsZero()) {
      if (image->Intrinsics() != intrinsics_matrix) {
        intrinsics_matrix = image->Intrinsics();
        fprintf(fp, "intrinsics_matrix  %g %g %g  %g %g %g  %g %g %g\n",
          intrinsics_matrix[0][0], intrinsics_matrix[0][1], intrinsics_matrix[0][2],
          intrinsics_matrix[1][0], intrinsics_matrix[1][1], intrinsics_matrix[1][2],
          intrinsics_matrix[2][0], intrinsics_matrix[2][1], intrinsics_matrix[2][2]);
      }
    }

    // Write images
    // if (!image->WriteChannels()) return 0;

    // Write image with extrinsics
    R4Matrix m = image->CameraToWorld().Matrix();
    const char *depth_filename = (image->DepthFilename()) ? image->DepthFilename() : "-";
    const char *color_filename = (image->ColorFilename()) ? image->ColorFilename() : "-";
    fprintf(fp, "image %s %s  %g %g %g %g  %g %g %g %g  %g %g %g %g  %g %g %g %g\n",
       depth_filename, color_filename, 
       m[0][0], m[0][1], m[0][2], m[0][3],
       m[1][0], m[1][1], m[1][2], m[1][3], 
       m[2][0], m[2][1], m[2][2], m[2][3],
       m[3][0], m[3][1], m[3][2], m[3][3]);
  }

  // Return success
  return 1;
}



////////////////////////////////////////////////////////////////////////
// OBJ output functions
////////////////////////////////////////////////////////////////////////

static int
WriteMtlFile(const RGBDConfiguration *configuration, const char *filename)
{
  // Check texture directory name
  if (!configuration->TextureDirectory()) return 1;

  // Open mtl file
  FILE *fp = fopen(filename, "w");
  if (!fp) {
    RNFail("Unable to open file %s", filename);
    return 0;
  }

  // Write materials
  for (int i = 0; i < configuration->NSurfaces(); i++) {
    RGBDSurface *surface = configuration->Surface(i);
    if (!surface->TextureFilename()) continue;
    fprintf(fp, "newmtl surface%d\n", i);
    fprintf(fp, "Kd 1 1 1\n");
    fprintf(fp, "map_Kd %s/%s\n", configuration->TextureDirectory(), surface->TextureFilename());
    fprintf(fp, "\n");
  }
  
  // Close mtl file
  fclose(fp);

  // Return success
  return 1;
}



int RGBDConfiguration::
WriteObjFile(const char *filename) const
{
  // Write mtl file
  char mtl_filename[1024];
  strncpy(mtl_filename, filename, 1020);
  char *endp = strrchr(mtl_filename, '.');
  if (endp) { *endp = '\0'; strcat(mtl_filename, ".mtl"); }
  if (!WriteMtlFile(this, mtl_filename)) return 0;

  // Open obj file
  FILE *fp = fopen(filename, "w");
  if (!fp) {
    RNFail("Unable to open file %s", filename);
    return 0;
  }

  // Write material library
  char *startp = strrchr(mtl_filename, '/');
  if (startp) startp++;
  else startp = mtl_filename;
  fprintf(fp, "mtllib %s\n", startp);

  // Write surfaces
  for (int s = 0; s < NSurfaces(); s++) {
    RGBDSurface *surface = Surface(s);
    int width = surface->NTexels(RN_X);
    int height = surface->NTexels(RN_Y);
    if ((width == 0) || (height == 0)) continue;
    
    // Write mesh
    R3Mesh *mesh = surface->mesh;
    if (mesh) {
      // Write material reference
      fprintf(fp, "g surface%d\n", s);
      fprintf(fp, "usemtl surface%d\n", s);

      // Write vertices
      for (int i = 0; i < mesh->NVertices(); i++) {
        R3MeshVertex *vertex = mesh->Vertex(i);
        const R3Point& p = mesh->VertexPosition(vertex);
        R2Point t, s = mesh->VertexTextureCoords(vertex);
        surface->TransformSurfaceToTexture(s, t);
        fprintf(fp, "vt %g %g\n", t.X() / width, t.Y() / height);
        fprintf(fp, "v %g %g %g\n", p.X(), p.Y(), p.Z());
      }

      // Write faces
      for (int i = 0; i < mesh->NFaces(); i++) {
        R3MeshFace *face = mesh->Face(i);
        R3MeshVertex *v0 = mesh->VertexOnFace(face, 0);
        R3MeshVertex *v1 = mesh->VertexOnFace(face, 1);
        R3MeshVertex *v2 = mesh->VertexOnFace(face, 2);
        int i0 = mesh->VertexID(v0) + 1;
        int i1 = mesh->VertexID(v1) + 1;
        int i2 = mesh->VertexID(v2) + 1;
        fprintf(fp, "f %d/%d %d/%d %d/%d\n", i0, i0, i1, i1, i2, i2);
      }
    }
  }
  
  // Close obj file
  fclose(fp);

  // Return success
  return 1;
}



////////////////////////////////////////////////////////////////////////
// Read/release functions
////////////////////////////////////////////////////////////////////////

int RGBDConfiguration::
ReadChannels(void)
{
  // Read images
  for (int i = 0; i < NImages(); i++) {
    RGBDImage *image = Image(i);
    image->ReadChannels();
  }

  // Read surfaces
  for (int i = 0; i < NSurfaces(); i++) {
    RGBDSurface *surface = Surface(i);
    surface->ReadChannels();
  }

  // Return success
  return 1;
}



int RGBDConfiguration::
ReleaseChannels(void)
{
  // Release images
  for (int i = 0; i < NImages(); i++) {
    RGBDImage *image = Image(i);
    image->ReleaseChannels();
  }

  // Release surfaces
  for (int i = 0; i < NSurfaces(); i++) {
    RGBDSurface *surface = Surface(i);
    surface->ReleaseChannels();
  }

  // Return success
  return 1;
}



int RGBDConfiguration::
ReadColorChannels(void)
{
  // Read images
  for (int i = 0; i < NImages(); i++) {
    RGBDImage *image = Image(i);
    image->ReadColorChannels();
  }

  // Read surfaces
  for (int i = 0; i < NSurfaces(); i++) {
    RGBDSurface *surface = Surface(i);
    surface->ReadColorChannels();
  }

  // Return success
  return 1;
}



int RGBDConfiguration::
ReleaseColorChannels(void)
{
  // Release images
  for (int i = 0; i < NImages(); i++) {
    RGBDImage *image = Image(i);
    image->ReleaseColorChannels();
  }

  // Release surfaces
  for (int i = 0; i < NSurfaces(); i++) {
    RGBDSurface *surface = Surface(i);
    surface->ReleaseColorChannels();
  }

  // Return success
  return 1;
}



int RGBDConfiguration::
ReadDepthChannels(void)
{
  // Read images
  for (int i = 0; i < NImages(); i++) {
    RGBDImage *image = Image(i);
    image->ReadDepthChannel();
  }

  // Return success
  return 1;
}



int RGBDConfiguration::
ReleaseDepthChannels(void)
{
  // Release images
  for (int i = 0; i < NImages(); i++) {
    RGBDImage *image = Image(i);
    image->ReleaseDepthChannel();
  }

  // Return success
  return 1;
}



int RGBDConfiguration::
ReadCategoryChannels(void)
{
  // Read images
  for (int i = 0; i < NImages(); i++) {
    RGBDImage *image = Image(i);
    image->ReadCategoryChannel();
  }

  // Return success
  return 1;
}



int RGBDConfiguration::
ReleaseCategoryChannels(void)
{
  // Release images
  for (int i = 0; i < NImages(); i++) {
    RGBDImage *image = Image(i);
    image->ReleaseCategoryChannel();
  }

  // Return success
  return 1;
}



int RGBDConfiguration::
ReadInstanceChannels(void)
{
  // Read images
  for (int i = 0; i < NImages(); i++) {
    RGBDImage *image = Image(i);
    image->ReadInstanceChannel();
  }

  // Return success
  return 1;
}



int RGBDConfiguration::
ReleaseInstanceChannels(void)
{
  // Release images
  for (int i = 0; i < NImages(); i++) {
    RGBDImage *image = Image(i);
    image->ReleaseInstanceChannel();
  }

  // Return success
  return 1;
}



////////////////////////////////////////////////////////////////////////
// Update functions
////////////////////////////////////////////////////////////////////////

void RGBDConfiguration::
InvalidateWorldBBox(void)
{
  // Mark bounding box for recomputation
  world_bbox.Reset(R3Point(FLT_MAX, FLT_MAX, FLT_MAX), R3Point(-FLT_MAX, -FLT_MAX, -FLT_MAX));
}



void RGBDConfiguration::
UpdateWorldBBox(void)
{
  // Initialize bounding box
  world_bbox = R3null_box;

  // Union image bounding boxes
  for (int i = 0; i < NImages(); i++) {
    RGBDImage *image = Image(i);
    world_bbox.Union(image->WorldViewpoint());
    if (!image->world_bbox.IsEmpty()) {
      world_bbox.Union(image->world_bbox);
    }
  }

  // Union surface bounding boxes
  for (int i = 0; i < NSurfaces(); i++) {
    RGBDSurface *surface = Surface(i);
    world_bbox.Union(surface->WorldBBox());
  }
}



RGBDImage *RGBDConfiguration::
AllocateImage(void) const
{
  // Allocate image (can be over-ridden by derived class)
  return new RGBDImage();
}



RGBDImage *RGBDConfiguration::
FindImage(const char *name) const
{
  // Find image with matching name
  for (int i = 0; i < NImages(); i++) {
    RGBDImage *image = Image(i);
    if (image->Name() && !strcmp(image->Name(), name)) return image;
  }

  // Image not found
  return NULL;
}



} // namespace gaps
