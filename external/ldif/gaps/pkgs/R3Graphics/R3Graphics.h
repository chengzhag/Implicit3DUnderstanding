/* Include file for R3 graphics module */
#ifndef __R3__GRAPHICS__H__
#define __R3__GRAPHICS__H__



/* Class declarations */

namespace gaps {
class R3Scene;
class R3SceneNode;
class R3SceneElement;
}



/* Dependency include files */

#include "R3Shapes/R3Shapes.h"



/* Material include files */

#include "R3Brdf.h"
#include "R2Texture.h"
#include "R3Material.h"



/* Light include files */

#include "R3Light.h"
#include "R3DirectionalLight.h"
#include "R3PointLight.h"
#include "R3SpotLight.h"
#include "R3AreaLight.h"



/* Viewing include files */

#include "R2Viewport.h"
#include "R3Camera.h"
#include "R3Frustum.h"
#include "R3Viewer.h"



/* Scene include files */

#include "R3SceneReference.h"
#include "R3SceneElement.h"
#include "R3SceneNode.h"
#include "R3Scene.h"



/* Initialization functions */

namespace gaps{
int R3InitGraphics(void);
void R3StopGraphics(void);
}



#endif


