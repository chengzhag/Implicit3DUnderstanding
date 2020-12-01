// Source file for the scalar image processing program



// Include files 

namespace gaps {}
using namespace gaps;
#include "R2Shapes/R2Shapes.h"
#include "RNMath/RNMath.h"



// Type definitions

typedef enum {
  NOP_OPERATION,
  ABS_OPERATION,
  SQUARE_OPERATION,
  SQRT_OPERATION,
  NEGATE_OPERATION,
  INVERT_OPERATION,
  NORMALIZE_OPERATION,
  GRADIENT_X_OPERATION,
  GRADIENT_Y_OPERATION,
  GRADIENT_MAGNITUDE_OPERATION,
  GRADIENT_ANGLE_OPERATION,
  HESSIAN_XX_OPERATION,
  HESSIAN_XY_OPERATION,
  HESSIAN_YX_OPERATION,
  HESSIAN_YY_OPERATION,
  LAPLACIAN_OPERATION,
  LAPLACIAN_X_OPERATION,
  LAPLACIAN_Y_OPERATION,
  HARRIS_CORNER_FILTER_OPERATION,
  DETECT_EDGES_OPERATION,
  DETECT_CORNERS_OPERATION,
  SIGNED_DISTANCE_OPERATION,
  SQUARED_DISTANCE_OPERATION,
  POINT_SYMMETRY_OPERATION,
  CLEAR_OPERATION,
  FILL_HOLES_OPERATION,
  FILL_GAPS_OPERATION,
  SUBSTITUTE_OPERATION,
  ADD_OPERATION,
  SUBTRACT_OPERATION,
  MULTIPLY_OPERATION,
  DIVIDE_OPERATION,
  POW_OPERATION,
  BILATERAL_FILTER_OPERATION,
  MIN_OPERATION,
  MAX_OPERATION,
  MEDIAN_OPERATION,
  PERCENTILE_OPERATION,
  MASK_NONMINIMA_OPERATION,
  MASK_NONMAXIMA_OPERATION,
  FLATTEN_OPERATION,
  CONNECTED_COMPONENT_OPERATION,
  CONNECTED_COMPONENT_SIZE_OPERATION,
  CONNECTED_COMPONENT_LABEL_OPERATION,
  DILATE_OPERATION,
  ERODE_OPERATION,
  BLUR_OPERATION,
  BLUR_X_OPERATION,
  BLUR_Y_OPERATION,
  ADD_NOISE_OPERATION,
  THRESHOLD_OPERATION,
  CROP_OPERATION,
  RESAMPLE_OPERATION,
  REFINE_OPERATION,
  ADD_GRID_OPERATION,
  SUBTRACT_GRID_OPERATION,
  MULTIPLY_GRID_OPERATION,
  DIVIDE_GRID_OPERATION,
  MASK_GRID_OPERATION,
  OVERLAY_GRID_OPERATION,
  THRESHOLD_GRID_OPERATION,
  SUN3D_BITSHIFT_OPERATION,
  POISSON_OPERATION,
  HASH_OPERATION,
  DU_OPERATION,
  DV_OPERATION,
  TMP_OPERATION,
  NUM_OPERATIONS
} OperationType;

struct Operation {
  int type;
  const char *operand1;
  const char *operand2;
  const char *operand3;
  const char *operand4;
  const char *operand5;
  const char *operand6;
};



// Program variables

static char *input_name = NULL;
static char *output_name = NULL;
static const int max_operations = 100;
static Operation operations[max_operations];
static int noperations = 0;
static int print_verbose = 0;
static int print_debug = 0;



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



static int 
WriteGrid(R2Grid *grid, const char *filename)
{
  // Start statistics
  RNTime start_time;
  start_time.Read();

  // Write grid
  int status = grid->Write(filename);

  // Print statistics
  if (print_verbose) {
    printf("Wrote grid to %s\n", filename);
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

  // Return status
  return status;
}



struct GridValue {
  int index;
  RNScalar delta;
  GridValue **heap_ptr;
};



RNScalar 
Delta(R2Grid *grid, int x, int y)
{
  // Check value
  RNScalar value = grid->GridValue(x, y);
  if (value == R2_GRID_UNKNOWN_VALUE) return R2_GRID_UNKNOWN_VALUE;

  // Find sum of neighbor values
  int neighbor_count = 0;
  RNScalar neighbor_sum = 0;
  RNScalar neighbor_minimum = FLT_MAX;
  for (int dy = -1; dy <= 1; dy++) {
    for (int dx = -1; dx <= 1; dx++) {
      if ((dx == 0) && (dy == 0)) continue;
      int i = x + dx;
      if ((i < 0) || (i >= grid->XResolution())) continue;
      int j = y + dy;
      if ((j < 0) || (j >= grid->YResolution())) continue;
      RNScalar neighbor_value = grid->GridValue(i, j);
      if (neighbor_value == R2_GRID_UNKNOWN_VALUE) continue;
      if (neighbor_value < neighbor_minimum) neighbor_minimum = neighbor_value;
      neighbor_sum += neighbor_value;
      neighbor_count++;
    }
  }

  // Check if has any neighbors
  if (neighbor_count == 0) return R2_GRID_UNKNOWN_VALUE;

  // Compute mean of neighbor values
  RNScalar neighbor_mean = neighbor_sum / neighbor_count;
  //RNScalar delta = value - neighbor_minimum;
  RNScalar delta = value - neighbor_mean;
  return delta;
}
    


void 
FlattenGrid(R2Grid *grid, RNScalar delta, RNScalar sigma)
{
  // Allocate array of grid value wrappers
  GridValue **gvs = new GridValue * [ grid->NEntries() ];
  for (int i = 0; i < grid->NEntries(); i++) gvs[i] = NULL;

  // Build priority queue of grid values
  int x, y;
  GridValue tmp;
  RNHeap<GridValue *> heap(&tmp, &tmp.delta, &tmp.heap_ptr);
  for (int index = 0; index < grid->NEntries(); index++) {
    RNScalar value = grid->GridValue(index);
    if (value != R2_GRID_UNKNOWN_VALUE) {
      grid->IndexToIndices(index, x, y);
      RNScalar target_value = grid->GridValue(x, y, sigma);
      if (target_value != R2_GRID_UNKNOWN_VALUE) {
        RNScalar d = target_value - value;
        GridValue *gv = new GridValue();
        gv->index = index;
        gv->delta = d;
        gv->heap_ptr = NULL;
        gvs[index] = gv;
        heap.Push(gv);
      }
    }
  }

  // Iteratively adjust grid values from largest negative step to smallest
  while (!heap.IsEmpty()) {
    // Get grid entry with highest delta
    GridValue *gv = heap.Peek();
    if (gv->delta > -delta) break;
    RNScalar value = grid->GridValue(gv->index);
    assert(value != R2_GRID_UNKNOWN_VALUE);

    // Update grid entry
    grid->SetGridValue(gv->index, value + gv->delta);

    // Update heap deltas for entry and all neighbors
    int x, y, index;
    grid->IndexToIndices(gv->index, x, y);
    // printf("%d %d %.6f\n", x, y, gv->delta);
    for (int dy = -1; dy <= 1; dy++) {
      for (int dx = -1; dx <= 1; dx++) {
        int i = x + dx;
        if ((i < 0) || (i >= grid->XResolution())) continue;
        int j = y + dy;
        if ((j < 0) || (j >= grid->YResolution())) continue;
        grid->IndicesToIndex(i, j, index);
        if (!gvs[index]) continue;
        RNScalar current_value = grid->GridValue(i, j);
        RNScalar target_value = grid->GridValue(i, j, sigma);
        assert(current_value != R2_GRID_UNKNOWN_VALUE);
        assert(target_value != R2_GRID_UNKNOWN_VALUE);
        gvs[index]->delta = target_value - current_value;
        heap.Update(gvs[index]);
      }
    }
  }

  // Delete memory
  for (int i = 0; i < grid->NEntries(); i++) { 
    if (gvs[i]) delete gvs[i];
  }
}


void 
Poisson(R2Grid *grid, R2Grid *dx_grid, R2Grid *dy_grid)
{
  // Create system of equations
  int n = grid->NEntries();
  int nx = grid->XResolution();
  RNSystemOfEquations equations(n);
  
  // Add inertia equations
  RNScalar inertia_weight = 1.0E-3;
  for (int i = 0; i < n; i++) {
    RNScalar d = grid->GridValue(i);
    RNPolynomial *f = new RNPolynomial(1.0, i, 1.0);
    RNAlgebraic *e = new RNAlgebraic(RN_SUBTRACT_OPERATION, f, d);
    e->Multiply(inertia_weight);
    equations.InsertEquation(e);
  }

  // Add dx equations
  RNScalar dx_weight = 1.0;
  for (int iy = 1; iy < grid->YResolution()-1; iy++) {
    for (int ix = 1; ix < grid->XResolution()-1; ix++) {
      RNScalar dx = dx_grid->GridValue(ix, iy);
      RNPolynomial *e = new RNPolynomial();
      e->AddTerm(-0.25, (iy-1)*nx+(ix-1), 1.0);
      e->AddTerm(-0.50, (iy)*nx+(ix-1), 1.0);
      e->AddTerm(-0.25, (iy+1)*nx+(ix-1), 1.0);
      e->AddTerm( 0.25, (iy-1)*nx+(ix+1), 1.0);
      e->AddTerm( 0.50, (iy)*nx+(ix+1), 1.0);
      e->AddTerm( 0.25, (iy+1)*nx+(ix+1), 1.0);
      e->Subtract(dx);
      e->Multiply(dx_weight);
      equations.InsertEquation(e);
    }
  }
  
  // Add dy equations
  RNScalar dy_weight = 1.0;
  for (int iy = 1; iy < grid->YResolution()-1; iy++) {
    for (int ix = 1; ix < grid->XResolution()-1; ix++) {
      RNScalar dy = dy_grid->GridValue(ix, iy);
      RNPolynomial *e = new RNPolynomial();
      e->AddTerm(-0.25, (iy-1)*nx+(ix-1), 1.0);
      e->AddTerm(-0.50, (iy-1)*nx+(ix), 1.0);
      e->AddTerm(-0.25, (iy-1)*nx+(ix+1), 1.0);
      e->AddTerm( 0.25, (iy+1)*nx+(ix-1), 1.0);
      e->AddTerm( 0.50, (iy+1)*nx+(ix), 1.0);
      e->AddTerm( 0.25, (iy+1)*nx+(ix+1), 1.0);
      e->Subtract(dy);
      e->Multiply(dy_weight);
      equations.InsertEquation(e);
    }
  }
  
  // Initialize variables
  double *x = new double [ n ];
  for (int i = 0; i < n; i++) x[i] = grid->GridValue(i);

  // Solve system of equations
  if (equations.NEquations() >= n) {
    if (!equations.Minimize(x, RN_CSPARSE_SOLVER, 1E-3)) {
      RNFail("Unable to minimize system of equations\n");
      delete [] x;
      return;
    }
  }

  // Copy solution into grid
  for (int i = 0; i < n; i++) grid->SetGridValue(i, x[i]);

  // Delete variables
  delete [] x;
}


struct LineSegment {
  LineSegment(R2Point p1, R2Point p2, RNScalar score) : p1(p1), p2(p2), score(score) {};
  R2Point p1, p2;
  RNScalar score;
};


#if 0
static int 
RNCompareScalarPtrs(const void *value1, const void *value2)
{
  const RNScalar **scalar1pp = (const RNScalar **) value1;
  const RNScalar **scalar2pp = (const RNScalar **) value2;
  const RNScalar *scalar1p = *scalar1pp;
  const RNScalar *scalar2p = *scalar2pp;
  if (*scalar1p < *scalar2p) return -1;
  else if (*scalar1p > *scalar2p) return 1;
  else return 0;
}
#endif



void SUN3DBitshift(R2Grid *grid)
{
  for (int i = 0; i < grid->NEntries(); i++) {
    unsigned int d = (unsigned int) (grid->GridValue(i) + 0.5);
    d = ((d >> 3) & 0x1FFF) | ((d & 0x7) << 13);
    grid->SetGridValue(i, d);
  }
}



void Hash(R2Grid *grid, int nvalues)
{
  // Get some stats
  RNInterval range = grid->Range();
  if (range.Diameter() == 0) return;
  RNScalar scale = 1.0 / range.Diameter();
  
  // Replace values with hash to value in [0:nvalues-1]
  for (int iy = 0; iy < grid->YResolution(); iy++) {
    for (int ix = 0; ix < grid->XResolution(); ix++) {
      RNScalar value = grid->GridValue(ix, iy);
      if (value == R2_GRID_UNKNOWN_VALUE) continue;
      int hash = (int) (nvalues * scale * (value - range.Min()));
      if (hash > nvalues-1) hash = nvalues-1;
      grid->SetGridValue(ix, iy, hash);
    }
  }
}



void DU(R2Grid *grid)
{
  R2Grid f(*grid);
  grid->Clear(R2_GRID_UNKNOWN_VALUE);
  for (int iy = 0; iy < grid->YResolution(); iy++) {
    for (int ix = 0; ix < grid->XResolution()-1; ix++) {
      RNScalar d0 = f.GridValue(ix, iy);
      if (RNIsZero(d0) || (d0 == R2_GRID_UNKNOWN_VALUE)) continue;
      RNScalar d1 = f.GridValue(ix+1, iy);
      if (RNIsZero(d1) || (d1 == R2_GRID_UNKNOWN_VALUE)) continue;
      grid->SetGridValue(ix, iy, d1 - d0);
    }
  }
}



void DV(R2Grid *grid)
{
  R2Grid f(*grid);
  grid->Clear(R2_GRID_UNKNOWN_VALUE);
  for (int iy = 0; iy < grid->YResolution()-1; iy++) {
    for (int ix = 0; ix < grid->XResolution(); ix++) {
      RNScalar d0 = f.GridValue(ix, iy);
      if (RNIsZero(d0) || (d0 == R2_GRID_UNKNOWN_VALUE)) continue;
      RNScalar d1 = f.GridValue(ix, iy+1);
      if (RNIsZero(d1) || (d1 == R2_GRID_UNKNOWN_VALUE)) continue;
      grid->SetGridValue(ix, iy, d1 - d0);
    }
  }
}



void Tmp(R2Grid *grid)
{
#if 0
  #include "fft/fft.h"
  R2Grid f(*grid);
  printf("%g ", f.L1Norm());
  R2ForwardFFT(f.XResolution(), f.YResolution(), (double *) f.GridValues());
  printf("%g ", f.L1Norm());
  f.Write("forward.pfm");
  f.Write("forward.ppm");
  R2ReverseFFT(f.XResolution(), f.YResolution(), (double *) f.GridValues());
  f.Write("reverse.pfm");
  printf("%g ", f.L1Norm());
  R2Grid diff(*grid);
  diff.Subtract(f);
  printf("%g\n", diff.L1Norm() / grid->L1Norm());
  *grid = f;
#endif
}



static int
ApplyOperations(R2Grid *grid, Operation *operations, int noperations)
{
  // Start statistics
  RNTime start_time;
  start_time.Read();

  // Apply operations
  for (int i = 0; i < noperations; i++) {
    Operation *operation = &operations[i];

    // Read grid operand
    R2Grid *grid1 = NULL;
    R2Grid *grid2 = NULL;
    switch (operation->type) {
    case ADD_GRID_OPERATION: 
    case SUBTRACT_GRID_OPERATION: 
    case MULTIPLY_GRID_OPERATION: 
    case DIVIDE_GRID_OPERATION: 
    case MASK_GRID_OPERATION: 
    case OVERLAY_GRID_OPERATION: 
    case THRESHOLD_GRID_OPERATION: 
      grid1 = ReadGrid(operation->operand1);
      if (!grid1) {
        RNFail("Unable to read grid file (%s) for operation %d\n", operation->operand1, i);
        return 0;
      }
      break;

    case POISSON_OPERATION: 
      grid1 = ReadGrid(operation->operand1);
      if (!grid1) {
        RNFail("Unable to read grid file (%s) for operation %d\n", operation->operand1, i);
        return 0;
      }
      grid2 = ReadGrid(operation->operand2);
      if (!grid2) {
        RNFail("Unable to read grid file (%s) for operation %d\n", operation->operand2, i);
        return 0;
      }
      break;
    }

    // Apply operation
    switch (operation->type) {
    case NOP_OPERATION: break;
    case ABS_OPERATION: grid->Abs(); break;
    case SQUARE_OPERATION: grid->Square(); break;
    case SQRT_OPERATION: grid->Sqrt(); break;
    case NEGATE_OPERATION: grid->Negate(); break;
    case INVERT_OPERATION: grid->Invert(); break;
    case NORMALIZE_OPERATION: grid->Normalize(); break;
    case DETECT_EDGES_OPERATION: grid->DetectEdges(); break;
    case DETECT_CORNERS_OPERATION: grid->DetectCorners(); break;
    case FILL_HOLES_OPERATION: grid->FillHoles(); break;
    case FILL_GAPS_OPERATION: grid->FillHoles(atoi(operation->operand1)); break;
    case GRADIENT_X_OPERATION: grid->Gradient(RN_X); break;
    case GRADIENT_Y_OPERATION: grid->Gradient(RN_Y); break;
    case GRADIENT_MAGNITUDE_OPERATION: grid->GradientMagnitude(); break;
    case GRADIENT_ANGLE_OPERATION: grid->GradientAngle(); break;
    case LAPLACIAN_OPERATION: grid->Laplacian(); break;
    case LAPLACIAN_X_OPERATION: grid->Laplacian(RN_X); break;
    case LAPLACIAN_Y_OPERATION: grid->Laplacian(RN_Y); break;
    case HARRIS_CORNER_FILTER_OPERATION: grid->HarrisCornerFilter(atoi(operation->operand1), atof(operation->operand2)); break;
    case HESSIAN_XX_OPERATION: grid->Hessian(RN_X, RN_X); break;
    case HESSIAN_XY_OPERATION: grid->Hessian(RN_X, RN_Y); break;
    case HESSIAN_YX_OPERATION: grid->Hessian(RN_Y, RN_X); break;
    case HESSIAN_YY_OPERATION: grid->Hessian(RN_Y, RN_Y); break;
    case SIGNED_DISTANCE_OPERATION: grid->SignedDistanceTransform(); break;
    case SQUARED_DISTANCE_OPERATION: grid->SquaredDistanceTransform(); break;
    case CLEAR_OPERATION: grid->Clear(atof(operation->operand1)); break;
    case ADD_OPERATION: grid->Add(atof(operation->operand1)); break;
    case SUBTRACT_OPERATION: grid->Subtract(atof(operation->operand1)); break;
    case MULTIPLY_OPERATION: grid->Multiply(atof(operation->operand1)); break;
    case DIVIDE_OPERATION: grid->Divide(atof(operation->operand1)); break;
    case POW_OPERATION: grid->Pow(atof(operation->operand1)); break;
    case BILATERAL_FILTER_OPERATION: grid->BilateralFilter(atof(operation->operand1), atof(operation->operand2)); break;
    case MAX_OPERATION: grid->MaxFilter(atof(operation->operand1)); break;
    case MIN_OPERATION: grid->MinFilter(atof(operation->operand1)); break;
    case MEDIAN_OPERATION: grid->MedianFilter(atof(operation->operand1)); break;
    case PERCENTILE_OPERATION: grid->PercentileFilter(atof(operation->operand1), atof(operation->operand2)); break;
    case MASK_NONMINIMA_OPERATION: grid->MaskNonMinima(atof(operation->operand1)); break;
    case MASK_NONMAXIMA_OPERATION: grid->MaskNonMaxima(atof(operation->operand1)); break;
    case FLATTEN_OPERATION: FlattenGrid(grid, atof(operation->operand1), atof(operation->operand2)); break;
    case CONNECTED_COMPONENT_SIZE_OPERATION: grid->ConnectedComponentSizeFilter(atof(operation->operand1)); break; 
    case CONNECTED_COMPONENT_LABEL_OPERATION: grid->ConnectedComponentLabelFilter(atof(operation->operand1)); break; 
    case DILATE_OPERATION: grid->Dilate(atof(operation->operand1)); break;
    case ERODE_OPERATION: grid->Erode(atof(operation->operand1)); break;
    case BLUR_OPERATION: grid->Blur(atof(operation->operand1)); break;
    case BLUR_X_OPERATION: grid->Blur(RN_X, atof(operation->operand1)); break;
    case BLUR_Y_OPERATION: grid->Blur(RN_Y, atof(operation->operand1)); break;
    case ADD_NOISE_OPERATION: grid->AddNoise(atof(operation->operand1)); break;
    case POINT_SYMMETRY_OPERATION: grid->PointSymmetryTransform(atoi(operation->operand1)); break;
    case CROP_OPERATION: *grid = R2Grid(*grid, atoi(operation->operand1), atoi(operation->operand2), atoi(operation->operand3), atoi(operation->operand4)); break;
    case RESAMPLE_OPERATION: grid->Resample(atoi(operation->operand1), atoi(operation->operand2)); break;
    case ADD_GRID_OPERATION: grid->Add(*grid1); break;
    case SUBTRACT_GRID_OPERATION: grid->Subtract(*grid1); break;
    case MULTIPLY_GRID_OPERATION: grid->Multiply(*grid1); break;
    case DIVIDE_GRID_OPERATION: grid->Divide(*grid1); break;
    case MASK_GRID_OPERATION: grid->Mask(*grid1); break;
    case OVERLAY_GRID_OPERATION: grid->Overlay(*grid1); break;
    case SUN3D_BITSHIFT_OPERATION: SUN3DBitshift(grid); break;
    case POISSON_OPERATION: Poisson(grid, grid1, grid2); break;
    case HASH_OPERATION: Hash(grid, atoi(operation->operand1)); break;
    case DU_OPERATION: DU(grid); break;
    case DV_OPERATION: DV(grid); break;
    case TMP_OPERATION: Tmp(grid); break;

    case REFINE_OPERATION: {
      RNScalar scale = atof(operation->operand1);
      int xres = (int) (scale * (grid->XResolution() - 1) + 1 + 0.5);
      int yres = (int) (scale * (grid->YResolution() - 1) + 1 + 0.5);
      grid->Resample(xres, yres);
      break; }

    case CONNECTED_COMPONENT_OPERATION: {
      RNScalar under_isolevel_value, too_small_value, too_large_value;
      if (!strcmp(operation->operand4, "keep")) under_isolevel_value = R2_GRID_KEEP_VALUE;
      else if (!strcmp(operation->operand4, "unknown")) under_isolevel_value = R2_GRID_UNKNOWN_VALUE;
      else under_isolevel_value = atof(operation->operand4);
      if (!strcmp(operation->operand5, "keep")) too_small_value = R2_GRID_KEEP_VALUE;
      else if (!strcmp(operation->operand5, "unknown")) too_small_value = R2_GRID_UNKNOWN_VALUE;
      else too_small_value = atof(operation->operand5);
      if (!strcmp(operation->operand6, "keep")) too_large_value = R2_GRID_KEEP_VALUE;
      else if (!strcmp(operation->operand5, "unknown")) too_large_value = R2_GRID_UNKNOWN_VALUE;
      else too_large_value = atof(operation->operand6);
      grid->ConnectedComponentFilter(atof(operation->operand1), atof(operation->operand2), atof(operation->operand3), under_isolevel_value, too_small_value, too_large_value);
      break; }

    case SUBSTITUTE_OPERATION: {
      RNScalar value1, value2;
      if (!strcmp(operation->operand1, "unknown")) value1 = R2_GRID_UNKNOWN_VALUE;
      else value1 = atof(operation->operand1);
      if (!strcmp(operation->operand2, "unknown")) value2 = R2_GRID_UNKNOWN_VALUE;
      else value2 = atof(operation->operand2);
      grid->Substitute(value1, value2); 
      break; }

    case THRESHOLD_OPERATION: 
    case THRESHOLD_GRID_OPERATION: {
      RNScalar value1, value2;
      if (!strcmp(operation->operand2, "keep")) value1 = R2_GRID_KEEP_VALUE;
      else if (!strcmp(operation->operand2, "input")) value1 = R2_GRID_INPUT_VALUE;
      else if (!strcmp(operation->operand2, "unknown")) value1 = R2_GRID_UNKNOWN_VALUE;
      else value1 = atof(operation->operand2);
      if (!strcmp(operation->operand3, "keep")) value2 = R2_GRID_KEEP_VALUE;
      else if (!strcmp(operation->operand3, "input")) value2 = R2_GRID_INPUT_VALUE;
      else if (!strcmp(operation->operand3, "unknown")) value2 = R2_GRID_UNKNOWN_VALUE;
      else value2 = atof(operation->operand3);
      if (operation->type == THRESHOLD_GRID_OPERATION) grid->Threshold(*grid1, value1, value2); 
      else grid->Threshold(atof(operation->operand1), value1, value2); 
      break; }

    default: 
      RNFail("Unknown operation type (%d) in operation %d\n", operation->type, i); 
      return 0; 
    }

    // Print debug message
    if (print_debug) {
      printf("Applied operation: %d %s %s %s %s %s %s\n", operation->type, 
             (operation->operand1) ? operation->operand1 : "-", 
             (operation->operand2) ? operation->operand2 : "-", 
             (operation->operand3) ? operation->operand3 : "-",
             (operation->operand4) ? operation->operand4 : "-",
             (operation->operand5) ? operation->operand5 : "-",
             (operation->operand6) ? operation->operand6 : "-");
      printf("  Cardinality = %d\n", grid->Cardinality());
      RNInterval grid_range = grid->Range();
      printf("  Minimum = %g\n", grid_range.Min());
      printf("  Maximum = %g\n", grid_range.Max());
      printf("  L1Norm = %g\n", grid->L1Norm());
      printf("  L2Norm = %g\n", grid->L2Norm());
      fflush(stdout);
    }
  }

  // Print statistics
  if (print_verbose) {
    printf("Applied operations ...\n");
    printf("  Time = %.2f seconds\n", start_time.Elapsed());
    printf("  # Operations = %d\n", noperations);
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
  return 1;
}



static int 
ParseArgs(int argc, char **argv)
{
  // Check number of arguments
  if (argc < 2) {
    printf("Usage: grd2grd inputfile outputfile [-v]\n");
    exit(0);
  }

  // Parse arguments
  argc--; argv++;
  while (argc > 0) {
    if ((*argv)[0] == '-') {
      if (!strcmp(*argv, "-v")) {
        print_verbose = 1; 
      }
      else if (!strcmp(*argv, "-debug")) {
        print_debug = 1; 
      }
      else if (!strcmp(*argv, "-abs")) {
        assert(noperations < max_operations);
        Operation *operation = &operations[noperations++];
        operation->type = ABS_OPERATION;
      }
      else if (!strcmp(*argv, "-square")) {
        assert(noperations < max_operations);
        Operation *operation = &operations[noperations++];
        operation->type = SQUARE_OPERATION;
      }
      else if (!strcmp(*argv, "-sqrt")) {
        assert(noperations < max_operations);
        Operation *operation = &operations[noperations++];
        operation->type = SQRT_OPERATION;
      }
      else if (!strcmp(*argv, "-negate")) {
        assert(noperations < max_operations);
        Operation *operation = &operations[noperations++];
        operation->type = NEGATE_OPERATION;
      }
      else if (!strcmp(*argv, "-invert")) {
        assert(noperations < max_operations);
        Operation *operation = &operations[noperations++];
        operation->type = INVERT_OPERATION;
      }
      else if (!strcmp(*argv, "-normalize")) {
        assert(noperations < max_operations);
        Operation *operation = &operations[noperations++];
        operation->type = NORMALIZE_OPERATION;
      }
      else if (!strcmp(*argv, "-edges")) {
        assert(noperations < max_operations);
        Operation *operation = &operations[noperations++];
        operation->type = DETECT_EDGES_OPERATION;
      }
      else if (!strcmp(*argv, "-corners")) {
        assert(noperations < max_operations);
        Operation *operation = &operations[noperations++];
        operation->type = DETECT_CORNERS_OPERATION;
      }
      else if (!strcmp(*argv, "-fill_holes")) {
        assert(noperations < max_operations);
        Operation *operation = &operations[noperations++];
        operation->type = FILL_HOLES_OPERATION;
      }
      else if (!strcmp(*argv, "-fill_gaps")) {
        assert(noperations < max_operations);
        Operation *operation = &operations[noperations++];
        operation->type = FILL_GAPS_OPERATION;
        argc--; argv++; operation->operand1 = *argv; 
      }
      else if (!strcmp(*argv, "-dx") || !strcmp(*argv, "-gradient_x")) {
        assert(noperations < max_operations);
        Operation *operation = &operations[noperations++];
        operation->type = GRADIENT_X_OPERATION;
      }
      else if (!strcmp(*argv, "-dy") || !strcmp(*argv, "-gradient_y")) {
        assert(noperations < max_operations);
        Operation *operation = &operations[noperations++];
        operation->type = GRADIENT_Y_OPERATION;
      }
      else if (!strcmp(*argv, "-gradient") || !strcmp(*argv, "-gradient_magnitude") || !strcmp(*argv, "-sobel")) {
        assert(noperations < max_operations);
        Operation *operation = &operations[noperations++];
        operation->type = GRADIENT_MAGNITUDE_OPERATION;
      }
      else if (!strcmp(*argv, "-gradient_angle") || !strcmp(*argv, "-gradient_direction")) {
        assert(noperations < max_operations);
        Operation *operation = &operations[noperations++];
        operation->type = GRADIENT_ANGLE_OPERATION;
      }
      else if (!strcmp(*argv, "-laplacian")) {
        assert(noperations < max_operations);
        Operation *operation = &operations[noperations++];
        operation->type = LAPLACIAN_OPERATION;
      }
      else if (!strcmp(*argv, "-lx") || !strcmp(*argv, "-laplacian_x")) {
        assert(noperations < max_operations);
        Operation *operation = &operations[noperations++];
        operation->type = LAPLACIAN_X_OPERATION;
      }
      else if (!strcmp(*argv, "-ly") || !strcmp(*argv, "-laplacian_y")) {
        assert(noperations < max_operations);
        Operation *operation = &operations[noperations++];
        operation->type = LAPLACIAN_Y_OPERATION;
      }
      else if (!strcmp(*argv, "-hessian_xx")) {
        assert(noperations < max_operations);
        Operation *operation = &operations[noperations++];
        operation->type = HESSIAN_XX_OPERATION;
      }
      else if (!strcmp(*argv, "-hessian_xy")) {
        assert(noperations < max_operations);
        Operation *operation = &operations[noperations++];
        operation->type = HESSIAN_XX_OPERATION;
      }
      else if (!strcmp(*argv, "-hessian_yx")) {
        assert(noperations < max_operations);
        Operation *operation = &operations[noperations++];
        operation->type = HESSIAN_YX_OPERATION;
      }
      else if (!strcmp(*argv, "-hessian_yy")) {
        assert(noperations < max_operations);
        Operation *operation = &operations[noperations++];
        operation->type = HESSIAN_YY_OPERATION;
      }
      else if (!strcmp(*argv, "-harris_corner_filter")) {
        assert(noperations < max_operations);
        Operation *operation = &operations[noperations++];
        argc--; argv++; operation->operand1 = *argv; 
        argc--; argv++; operation->operand2 = *argv; 
        operation->type = HARRIS_CORNER_FILTER_OPERATION;
      }
      else if (!strcmp(*argv, "-signed_distance")) {
        assert(noperations < max_operations);
        Operation *operation = &operations[noperations++];
        operation->type = SIGNED_DISTANCE_OPERATION;
      }
      else if (!strcmp(*argv, "-squared_distance")) {
        assert(noperations < max_operations);
        Operation *operation = &operations[noperations++];
        operation->type = SQUARED_DISTANCE_OPERATION;
      }
      else if (!strcmp(*argv, "-clear")) {
        assert(noperations < max_operations);
        Operation *operation = &operations[noperations++];
        operation->type = CLEAR_OPERATION;
        argc--; argv++; operation->operand1 = *argv; 
      }
      else if (!strcmp(*argv, "-substitute")) {
        assert(noperations < max_operations);
        Operation *operation = &operations[noperations++];
        operation->type = SUBSTITUTE_OPERATION;
        argc--; argv++; operation->operand1 = *argv; 
        argc--; argv++; operation->operand2 = *argv; 
      }
      else if (!strcmp(*argv, "-add")) {
        assert(noperations < max_operations);
        Operation *operation = &operations[noperations++];
        operation->type = ADD_OPERATION;
        argc--; argv++; operation->operand1 = *argv; 
      }
      else if (!strcmp(*argv, "-subtract")) {
        assert(noperations < max_operations);
        Operation *operation = &operations[noperations++];
        operation->type = SUBTRACT_OPERATION;
        argc--; argv++; operation->operand1 = *argv; 
      }
      else if (!strcmp(*argv, "-multiply")) {
        assert(noperations < max_operations);
        Operation *operation = &operations[noperations++];
        operation->type = MULTIPLY_OPERATION;
        argc--; argv++; operation->operand1 = *argv; 
      }
      else if (!strcmp(*argv, "-divide")) {
        assert(noperations < max_operations);
        Operation *operation = &operations[noperations++];
        operation->type = DIVIDE_OPERATION;
        argc--; argv++; operation->operand1 = *argv; 
      }
      else if (!strcmp(*argv, "-pow")) {
        assert(noperations < max_operations);
        Operation *operation = &operations[noperations++];
        operation->type = POW_OPERATION;
        argc--; argv++; operation->operand1 = *argv; 
      }
      else if (!strcmp(*argv, "-point_symmetry")) {
        assert(noperations < max_operations);
        Operation *operation = &operations[noperations++];
        operation->type = POINT_SYMMETRY_OPERATION;
        argc--; argv++; operation->operand1 = *argv; 
      }
      else if (!strcmp(*argv, "-bilateral")) {
        assert(noperations < max_operations);
        Operation *operation = &operations[noperations++];
        operation->type = BILATERAL_FILTER_OPERATION;
        argc--; argv++; operation->operand1 = *argv; 
        argc--; argv++; operation->operand2 = *argv; 
      }
      else if (!strcmp(*argv, "-min")) {
        assert(noperations < max_operations);
        Operation *operation = &operations[noperations++];
        operation->type = MIN_OPERATION;
        argc--; argv++; operation->operand1 = *argv; 
      }
      else if (!strcmp(*argv, "-max")) {
        assert(noperations < max_operations);
        Operation *operation = &operations[noperations++];
        operation->type = MAX_OPERATION;
        argc--; argv++; operation->operand1 = *argv; 
      }
      else if (!strcmp(*argv, "-median")) {
        assert(noperations < max_operations);
        Operation *operation = &operations[noperations++];
        operation->type = MEDIAN_OPERATION;
        argc--; argv++; operation->operand1 = *argv; 
      }
      else if (!strcmp(*argv, "-percentile")) {
        assert(noperations < max_operations);
        Operation *operation = &operations[noperations++];
        operation->type = PERCENTILE_OPERATION;
        argc--; argv++; operation->operand1 = *argv; 
        argc--; argv++; operation->operand2 = *argv; 
      }
      else if (!strcmp(*argv, "-mask_nonminima")) {
        assert(noperations < max_operations);
        Operation *operation = &operations[noperations++];
        operation->type = MASK_NONMINIMA_OPERATION;
        argc--; argv++; operation->operand1 = *argv; 
      }
      else if (!strcmp(*argv, "-mask_nonmaxima")) {
        assert(noperations < max_operations);
        Operation *operation = &operations[noperations++];
        operation->type = MASK_NONMAXIMA_OPERATION;
        argc--; argv++; operation->operand1 = *argv; 
      }
      else if (!strcmp(*argv, "-flatten")) {
        assert(noperations < max_operations);
        Operation *operation = &operations[noperations++];
        operation->type = FLATTEN_OPERATION;
        argc--; argv++; operation->operand1 = *argv; 
        argc--; argv++; operation->operand2 = *argv; 
      }
      else if (!strcmp(*argv, "-connected_component")) {
        assert(noperations < max_operations);
        Operation *operation = &operations[noperations++];
        operation->type = CONNECTED_COMPONENT_OPERATION;
        argc--; argv++; operation->operand1 = *argv; 
        argc--; argv++; operation->operand2 = *argv; 
        argc--; argv++; operation->operand3 = *argv; 
        argc--; argv++; operation->operand4 = *argv; 
        argc--; argv++; operation->operand5 = *argv; 
        argc--; argv++; operation->operand6 = *argv; 
      }
      else if (!strcmp(*argv, "-connected_component_size")) {
        assert(noperations < max_operations);
        Operation *operation = &operations[noperations++];
        operation->type = CONNECTED_COMPONENT_SIZE_OPERATION;
        argc--; argv++; operation->operand1 = *argv; 
      }
      else if (!strcmp(*argv, "-connected_component_label")) {
        assert(noperations < max_operations);
        Operation *operation = &operations[noperations++];
        operation->type = CONNECTED_COMPONENT_LABEL_OPERATION;
        argc--; argv++; operation->operand1 = *argv; 
      }
      else if (!strcmp(*argv, "-dilate")) {
        assert(noperations < max_operations);
        Operation *operation = &operations[noperations++];
        operation->type = DILATE_OPERATION;
        argc--; argv++; operation->operand1 = *argv; 
      }
      else if (!strcmp(*argv, "-erode")) {
        assert(noperations < max_operations);
        Operation *operation = &operations[noperations++];
        operation->type = ERODE_OPERATION;
        argc--; argv++; operation->operand1 = *argv; 
      }
      else if (!strcmp(*argv, "-blur")) {
        assert(noperations < max_operations);
        Operation *operation = &operations[noperations++];
        operation->type = BLUR_OPERATION;
        argc--; argv++; operation->operand1 = *argv; 
      }
      else if (!strcmp(*argv, "-blur_x")) {
        assert(noperations < max_operations);
        Operation *operation = &operations[noperations++];
        operation->type = BLUR_X_OPERATION;
        argc--; argv++; operation->operand1 = *argv; 
      }
      else if (!strcmp(*argv, "-blur_y")) {
        assert(noperations < max_operations);
        Operation *operation = &operations[noperations++];
        operation->type = BLUR_Y_OPERATION;
        argc--; argv++; operation->operand1 = *argv; 
      }
      else if (!strcmp(*argv, "-add_noise")) {
        assert(noperations < max_operations);
        Operation *operation = &operations[noperations++];
        operation->type = ADD_NOISE_OPERATION;
        argc--; argv++; operation->operand1 = *argv; 
      }
      else if (!strcmp(*argv, "-threshold")) {
        assert(noperations < max_operations);
        Operation *operation = &operations[noperations++];
        operation->type = THRESHOLD_OPERATION;
        argc--; argv++; operation->operand1 = *argv; 
        argc--; argv++; operation->operand2 = *argv; 
        argc--; argv++; operation->operand3 = *argv; 
      }
      else if (!strcmp(*argv, "-crop")) {
        assert(noperations < max_operations);
        Operation *operation = &operations[noperations++];
        operation->type = CROP_OPERATION;
        argc--; argv++; operation->operand1 = *argv; 
        argc--; argv++; operation->operand2 = *argv; 
        argc--; argv++; operation->operand3 = *argv; 
        argc--; argv++; operation->operand4 = *argv; 
      }
      else if (!strcmp(*argv, "-resample")) {
        assert(noperations < max_operations);
        Operation *operation = &operations[noperations++];
        operation->type = RESAMPLE_OPERATION;
        argc--; argv++; operation->operand1 = *argv; 
        argc--; argv++; operation->operand2 = *argv; 
      }
      else if (!strcmp(*argv, "-refine")) {
        assert(noperations < max_operations);
        Operation *operation = &operations[noperations++];
        operation->type = REFINE_OPERATION;
        argc--; argv++; operation->operand1 = *argv; 
      }
      else if (!strcmp(*argv, "-add_grid")) {
        assert(noperations < max_operations);
        Operation *operation = &operations[noperations++];
        operation->type = ADD_GRID_OPERATION;
        argc--; argv++; operation->operand1 = *argv;
      }
      else if (!strcmp(*argv, "-subtract_grid")) {
        assert(noperations < max_operations);
        Operation *operation = &operations[noperations++];
        operation->type = SUBTRACT_GRID_OPERATION;
        argc--; argv++; operation->operand1 = *argv;
      }
      else if (!strcmp(*argv, "-multiply_grid")) {
        assert(noperations < max_operations);
        Operation *operation = &operations[noperations++];
        operation->type = MULTIPLY_GRID_OPERATION;
        argc--; argv++; operation->operand1 = *argv;
      }
      else if (!strcmp(*argv, "-divide_grid")) {
        assert(noperations < max_operations);
        Operation *operation = &operations[noperations++];
        operation->type = DIVIDE_GRID_OPERATION;
        argc--; argv++; operation->operand1 = *argv;
      }
      else if (!strcmp(*argv, "-mask_grid")) {
        assert(noperations < max_operations);
        Operation *operation = &operations[noperations++];
        operation->type = MASK_GRID_OPERATION;
        argc--; argv++; operation->operand1 = *argv;
      }
      else if (!strcmp(*argv, "-overlay_grid")) {
        assert(noperations < max_operations);
        Operation *operation = &operations[noperations++];
        operation->type = OVERLAY_GRID_OPERATION;
        argc--; argv++; operation->operand1 = *argv;
      }
      else if (!strcmp(*argv, "-threshold_grid")) {
        assert(noperations < max_operations);
        Operation *operation = &operations[noperations++];
        operation->type = THRESHOLD_GRID_OPERATION;
        argc--; argv++; operation->operand1 = *argv; 
        argc--; argv++; operation->operand2 = *argv; 
        argc--; argv++; operation->operand3 = *argv; 
      }
      else if (!strcmp(*argv, "-sun3d_bitshift")) {
        assert(noperations < max_operations);
        Operation *operation = &operations[noperations++];
        operation->type = SUN3D_BITSHIFT_OPERATION;
      }
      else if (!strcmp(*argv, "-poisson")) {
        assert(noperations < max_operations);
        Operation *operation = &operations[noperations++];
        operation->type = POISSON_OPERATION;
        argc--; argv++; operation->operand1 = *argv; 
        argc--; argv++; operation->operand2 = *argv; 
      }
      else if (!strcmp(*argv, "-hash")) {
        assert(noperations < max_operations);
        Operation *operation = &operations[noperations++];
        operation->type = HASH_OPERATION;
        argc--; argv++; operation->operand1 = *argv; 
      }
      else if (!strcmp(*argv, "-du")) {
        assert(noperations < max_operations);
        Operation *operation = &operations[noperations++];
        operation->type = DU_OPERATION;
      }
      else if (!strcmp(*argv, "-dv")) {
        assert(noperations < max_operations);
        Operation *operation = &operations[noperations++];
        operation->type = DV_OPERATION;
      }
      else if (!strcmp(*argv, "-tmp")) {
        assert(noperations < max_operations);
        Operation *operation = &operations[noperations++];
        operation->type = TMP_OPERATION;
      }
      else { 
        RNFail("Invalid program argument: %s", *argv); 
        exit(1); 
      }
    }
    else {
      if (!input_name) input_name = *argv;
      else if (!output_name) output_name = *argv;
      else { 
        RNFail("Invalid program argument: %s", *argv); 
        exit(1); 
      }
    }
    argv++; argc--;
  }

  // Check input filename
  if (!input_name) {
    RNFail("You did not specify an input file.\n");
    return 0;
  }

  // Check output filename
  if (!output_name) {
    RNFail("You did not specify an output file.\n");
    return 0;
  }

  // Return OK status 
  return 1;
}



int 
main(int argc, char **argv)
{
  // Parse program arguments
  if (!ParseArgs(argc, argv)) exit(-1);

  // Read grid file
  R2Grid *grid = ReadGrid(input_name);
  if (!grid) exit(-1);

  // Apply operations
  int status1 = ApplyOperations(grid, operations, noperations);
  if (!status1) exit(-1);

  // Write grid file
  int status2 = WriteGrid(grid, output_name);
  if (!status2) exit(-1);

  // Return success
  return 0;
}
