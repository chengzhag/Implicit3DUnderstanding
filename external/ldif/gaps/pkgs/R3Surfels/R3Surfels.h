/* Include file for R3 surfels module */

#ifndef __R3__SURFELS__H__
#define __R3__SURFELS__H__



/* Dependency include files */

#include "R3Graphics/R3Graphics.h"



/* Draw flags */

/* Shading draw flags */
#define R3_SURFEL_COLOR_DRAW_FLAG      0x0001
#define R3_SURFEL_NORMAL_DRAW_FLAG     0x0002
#define R3_SURFEL_IDENTIFIER_DRAW_FLAG 0x0004
/* Shape draw flags */
#define R3_SURFEL_DISC_DRAW_FLAG       0x0010
#define R3_SURFEL_DEFAULT_DRAW_FLAGS   R3_SURFEL_COLOR_DRAW_FLAG



/* Drawing method */

// Define only one of these
// #define R3_SURFEL_DRAW_WITH_DISPLAY_LIST
// #define R3_SURFEL_DRAW_WITH_VBO
// #define R3_SURFEL_DRAW_WITH_ARRAYS
#define R3_SURFEL_DRAW_WITH_GLBEGIN



/* Class declarations */

namespace gaps {
class R3Surfel;
class R3SurfelBlock;
class R3SurfelDatabase;
class R3SurfelConstraint;
class R3SurfelPoint;
class R3SurfelPointSet;
class R3SurfelPointGraph;
class R3SurfelNode;
class R3SurfelNodeSet;
class R3SurfelTree;
class R3SurfelScan;
class R3SurfelImage;
class R3SurfelFeature;
class R3SurfelFeatureSet;
class R3SurfelFeatureVector;
class R3SurfelObject;
class R3SurfelObjectSet;
class R3SurfelObjectProperty;
class R3SurfelObjectRelationship;
class R3SurfelLabel;
class R3SurfelLabelSet;
class R3SurfelLabelProperty;
class R3SurfelLabelRelationship;
class R3SurfelLabelAssignment;
typedef R3SurfelLabelAssignment R3SurfelObjectAssignment;
class R3SurfelScene;
}



/* Surfel pkg include files */

#include "R3Surfel.h"
#include "R3SurfelBlock.h"
#include "R3SurfelDatabase.h"
#include "R3SurfelConstraint.h"
#include "R3SurfelPoint.h"
#include "R3SurfelPointSet.h"
#include "R3SurfelPointGraph.h"
#include "R3SurfelNode.h"
#include "R3SurfelNodeSet.h"
#include "R3SurfelTree.h"
#include "R3SurfelScan.h"
#include "R3SurfelImage.h"
#include "R3SurfelFeature.h"
#include "R3SurfelFeatureSet.h"
#include "R3SurfelFeatureVector.h"
#include "R3SurfelFeatureEvaluation.h"
#include "R3SurfelObject.h"
#include "R3SurfelObjectSet.h"
#include "R3SurfelObjectProperty.h"
#include "R3SurfelObjectRelationship.h"
#include "R3SurfelLabel.h"
#include "R3SurfelLabelSet.h"
#include "R3SurfelLabelProperty.h"
#include "R3SurfelLabelRelationship.h"
#include "R3SurfelLabelAssignment.h"
#include "R3SurfelScene.h"



/* Utility include files */

#include "R3SurfelUtils.h"
#include "R3SurfelViewer.h"



/* Initialization functions */

namespace gaps {
int R3InitSurfels(void);
void R3StopSurfels(void);
}



#endif








