#include <iostream>
#include <glm/glm.hpp>
#include <SDL.h>
#include "SDLauxiliary.h"
#include "TestModelH.h"
#include <stdint.h>
#include <limits.h>

using namespace std;
using glm::vec3;
using glm::mat3;
using glm::vec4;
using glm::mat4;

SDL_Event event;

#define SCREEN_WIDTH 320
#define SCREEN_HEIGHT 256
#define FULLSCREEN_MODE false



/* ----------------------------------------------------------------------------*/
/* VARIABLES                                                                   */
/* ----------------------------------------------------------------------------*/

float focalLength = 100;

/* ----------------------------------------------------------------------------*/
/* STRUCTS                                                                   */
/* ----------------------------------------------------------------------------*/

struct Intersection {
       vec4 position;
       float distance;
       int triangleIndex;
};

/* ----------------------------------------------------------------------------*/
/* FUNCTIONS                                                                   */
/* ----------------------------------------------------------------------------*/


bool Update();
void Draw(screen* screen);
bool ClosestIntersection(
       vec4 start,
       vec4 dir,
       const vector<Triangle>& triangles,
       Intersection& closestIntersection );

int main( int argc, char* argv[] )
{

  screen *screen = InitializeSDL( SCREEN_WIDTH, SCREEN_HEIGHT, FULLSCREEN_MODE );


  while ( Update())
    {
      Draw(screen);
      SDL_Renderframe(screen);
    }

  SDL_SaveImage( screen, "screenshot.bmp" );

  KillSDL(screen);
  return 0;
}

/*Place your drawing here*/
void Draw(screen* screen) {
  /* Clear buffer */
  memset(screen->buffer, 0, screen->height*screen->width*sizeof(uint32_t));

  vec3 colour(1.0,0.0,0.0);
  for(int i=0; i<1000; i++)
    {
      uint32_t x = rand() % screen->width;
      uint32_t y = rand() % screen->height;
      PutPixelSDL(screen, x, y, colour);
    }

  // Raytrace begin
  vector<Triangle> triangles;
  LoadTestModel( triangles );

  for (size_t y = 0; y < SCREEN_HEIGHT; y++) {
    for (size_t x = 0; x < SCREEN_WIDTH; x++) {
      vec4 cameraPos( 100, 100, 20, 1.0);
      vec4 dir = vec4(x - (SCREEN_WIDTH/2), y - (SCREEN_HEIGHT/2), focalLength, 1.0);
      Intersection closestIntersection;

      if (ClosestIntersection(cameraPos, dir, triangles, closestIntersection) == true) {
        PutPixelSDL(screen, x, y, triangles[closestIntersection.triangleIndex].color);
      }
      else {
        PutPixelSDL(screen, x, y, colour);
      }
    }
  }
}

/*Place updates of parameters here*/
bool Update() {
  static int t = SDL_GetTicks();
  /* Compute frame time */
  int t2 = SDL_GetTicks();
  float dt = float(t2-t);
  t = t2;

  SDL_Event e;
  while(SDL_PollEvent(&e))
    {
      if (e.type == SDL_QUIT)
	{
	  return false;
	}
      else
	if (e.type == SDL_KEYDOWN)
	  {
	    int key_code = e.key.keysym.sym;
	    switch(key_code)
	      {
	      case SDLK_UP:
		/* Move camera forward */
		break;
	      case SDLK_DOWN:
		/* Move camera backwards */
		break;
	      case SDLK_LEFT:
		/* Move camera left */
		break;
	      case SDLK_RIGHT:
		/* Move camera right */
		break;
	      case SDLK_ESCAPE:
		/* Move camera quit */
		return false;
	      }
	  }
    }
  return true;
}

bool ClosestIntersection( vec4 start, vec4 dir,
                          const vector<Triangle>& triangles,
                          Intersection& closestIntersection ) {

  for (size_t i = 0; i < triangles.size(); i++) {
    vec4 v0 = triangles[i].v0;
    vec4 v1 = triangles[i].v1;
    vec4 v2 = triangles[i].v2;

    vec3 e1 = vec3(v1.x-v0.x,v1.y-v0.y,v1.z-v0.z);
    vec3 e2 = vec3(v2.x-v0.x,v2.y-v0.y,v2.z-v0.z);
    vec3 b = vec3(start.x-v0.x,start.y-v0.y,start.z-v0.z);

    // bottleneck
    // glm::mat3 A;
    vec3 dir3 = vec3(dir[0], dir[1], dir[2]);
    glm::mat3 A( -dir3, e1, e2 );

    vec3 x = glm::inverse( A ) * b;

    vec4 x4 = vec4(x[0], x[1], x[2], 1.0);

    // Don't check for intersection for large t
    float bound = std::numeric_limits<float>::max();

    // Check against inequalities, if true intersection occurred
    if (x[1] > 0 && x[2] > 0 && (x[1] + x[2]) < 1 and x[0] >= 0 and x[0] < bound) {
      closestIntersection.position = x4;
      closestIntersection.distance = x[0];
      closestIntersection.triangleIndex = i;
      return true;
    }
    return false;
  }
}
