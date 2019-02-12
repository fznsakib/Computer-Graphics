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

float focalLength = 20;

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

  vec3 colour(0.0,0.0,0.0);
  for(int i=0; i<1000; i++)
    {
      uint32_t x = rand() % screen->width;
      uint32_t y = rand() % screen->height;
      PutPixelSDL(screen, x, y, colour);
    }

  // Raytrace begin
  vector<Triangle> triangles;
  LoadTestModel( triangles );
  vec4 cameraPos( 0, 0, focalLength, 1.0);

  for (size_t v = -SCREEN_HEIGHT/2; v < SCREEN_HEIGHT/2; v++) {
    for (size_t u = -SCREEN_WIDTH/2; u < SCREEN_WIDTH/2; u++) {
      vec4 dir = vec4(u - (SCREEN_WIDTH/2), v - (SCREEN_HEIGHT/2), focalLength, 1.0);

      // vec4 normDir = 1/(sqrt(dir * dir)) * dir;
      Intersection closestIntersection;

      if (ClosestIntersection(cameraPos, dir, triangles, closestIntersection) == true) {
        // std::cout << "Intersection found, distance = " + std::to_string(closestIntersection.distance) + "\n";
        PutPixelSDL(screen, u, v, triangles[closestIntersection.triangleIndex].color);
      }
      else {
        // std::cout << "No intersection: paint black\n";
        PutPixelSDL(screen, u, v, colour);
      }

      // std::cout << std::to_string(x) + ", " + std::to_string(y) + " printed " << '\n';
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

    // Don't check for intersection for large t
    float bound = std::numeric_limits<float>::max();

    closestIntersection.distance = bound;

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

      bool check1 = x[1] > 0;
      bool check2 = x[2] > 0;
      bool check3 = (x[1] + x[2]) < 1;
      bool check4 = x[0] >= 0;
      bool check5 = x[0] < bound;
      bool check6 = x[0] < closestIntersection.distance;

      // std::cout << to_string(check1) + to_string(check2) + to_string(check3) + to_string(check4) + to_string(check5) + to_string(check6) << '\n';
      // Check against inequalities, if true intersection occurred
      if ( check1 && check2 && check3 && check4 && check5 && check6 ) {
        closestIntersection.position = x4;
        closestIntersection.distance = x[0];
        closestIntersection.triangleIndex = i;
        std::cout << "accept" << '\n';
        return true;
      }
    // else {
    //   return false;
    // }
      if ((i == triangles.size() - 1) && (closestIntersection.distance == bound)) {
        std::cout << "reject" << '\n';
        return false;
      }
      return false;
  }
}
