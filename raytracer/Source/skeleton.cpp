#include <iostream>
#include <glm/glm.hpp>
#include <SDL.h>
#include "SDLauxiliary.h"
#include "TestModelH.h"
#include <stdint.h>
#include <limits.h>
#include <math.h>

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

float focalLength = 128;
vec4 cameraPos( 0.0, 0.0, -2.0, 1.0);
float yaw = 0.0;
mat4 R(1.0f);

/* ----------------------------------------------------------------------------*/
/* STRUCTS                                                                   */
/* ----------------------------------------------------------------------------*/

struct Intersection {
       vec4 position;
       float distance;
       int triangleIndex;
};

struct Light {
      vec4 position;
      vec3 colour;
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
vec3 DirectLight( const Intersection &i, const vector<Triangle>& triangles, Light light );


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

  // Initialise triangles
  vector<Triangle> triangles;
  LoadTestModel( triangles );

  // Initialise lights
  vector<Light> lights;
  Light light1;
  light1.position = vec4( 0.0, -0.5, -0.7, 1.0 );
  light1.colour = vec3( 14.f * vec3( 1, 1, 1 ));
  lights.push_back(light1);

  // u and v are coordinates on the 2D screen
  for (int v = 0; v < SCREEN_HEIGHT; v++) {
    for (int u = 0; u < SCREEN_WIDTH; u++) {

      vec4 dir = (vec4((u - SCREEN_WIDTH/2), v - SCREEN_HEIGHT/2, focalLength, 1.0));
      /* CALCULATE THE NEW DIRECTION VECTOR AFTER YAW ROTATION */
      dir = R * dir;

      Intersection intersection;
      bool closestIntersection = ClosestIntersection(cameraPos, dir, triangles, intersection);

      if (closestIntersection == true) {
        // PutPixelSDL(screen, u, v, triangles[intersection.triangleIndex].color);

        vec3 pixelColour = vec3( 1.0, 1.0, 1.0 );
        for (int i = 0; i < lights.size(); i++) {
          pixelColour = DirectLight( intersection, triangles, lights[i] );
        }
        PutPixelSDL(screen, u, v, pixelColour);
      }
      else {
        PutPixelSDL(screen, u, v, colour);
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

  cout << "Render time : " << dt << " ms." <<endl;
  cout << "X : " + to_string(cameraPos[0]) + " Y : " + to_string(cameraPos[1]) + " Z : " + to_string(cameraPos[2]) + " Yaw : " + to_string(yaw) << endl;

  SDL_Event e;
  while(SDL_PollEvent(&e))
    {
      if (e.type == SDL_QUIT)
	    {
	      return false;
	    }
      else if (e.type == SDL_KEYDOWN)
	    {
	      int key_code = e.key.keysym.sym;
	      switch(key_code)
	      {
	        case SDLK_UP:
		        cameraPos += vec4(0, 0, 0.5, 0);
            cout << "KEYPRESS UP." <<endl;
		      break;
	        case SDLK_DOWN:
		        cameraPos += vec4(0, 0, -0.5, 0);
            cout << "KEYPRESS DOWN." <<endl;
		      break;
	        case SDLK_LEFT:
		        cameraPos += vec4(-0.5, 0, 0, 0);
            cout << "KEYPRESS <-." <<endl;
		      break;
	        case SDLK_RIGHT:
		        cameraPos += vec4(0.5, 0, 0, 0);
            cout << "KEYPRESS ->." <<endl;
		      break;
          case SDLK_a:
            cout << "KEYPRESS a." <<endl;
            yaw -= 0.174533;
            R[0][0] = cos(yaw);   R[0][1] = 0;   R[0][2] = -sin(yaw);
            R[1][0] = 0;          R[1][1] = 1;   R[1][2] = 0;
            R[2][0] = sin(yaw);   R[2][1] = 0;   R[2][2] = cos(yaw);
          break;

            /* POTENTIALLY NOT NEEDED ANYMORE */
            // rightAxis   = normalize(vec4(R[0][0], R[0][1], R[0][2], 1);
            // downAxis    = normalize(vec4(R[1][0], R[1][1], R[1][2], 1);
            // forwardAxis = normalize(vec4(R[2][0], R[2][1], R[2][2], 1);
            // vec4 right(cos(yaw), 0, -sin(yaw), 1);
            // vec4 down(0, 1, 0, 1);
            // vec4 forward(sin(yaw), 0, cos(yaw), 1);

          case SDLK_d:
            yaw += 0.174533;
            R[0][0] = cos(yaw);   R[0][1] = 0;   R[0][2] = -sin(yaw);
            R[1][0] = 0;          R[1][1] = 1;   R[1][2] = 0;
            R[2][0] = sin(yaw);   R[2][1] = 0;   R[2][2] = cos(yaw);
          break;

          case SDLK_i:
            focalLength += 10;
          break;

          case SDLK_o:
            focalLength -= 10;
          break;

	        case SDLK_ESCAPE:
            return false;
	      }
	    }
    }
  return true;
}


vec3 DirectLight( const Intersection &i, const vector<Triangle>& triangles, Light light ) {

  vec3 power;
  Triangle triangle = triangles[i.triangleIndex];
  vec4 r = sqrt((light.position - i.position) * (light.position - i.position));

  vec4 direction = light.position - i.position;
  direction = glm::normalize(direction);

  float a = glm::dot(direction, triangle.normal);
  float b = 4 * M_PI;

  vec3 surfaceArea = vec3(b * (r * r));

  // If light does not hit triangle
  if (a < 0) a = 0;

  power = (light.colour * a)/surfaceArea;

  return power;
}



bool ClosestIntersection( vec4 start, vec4 dir,
                          const vector<Triangle>& triangles,
                          Intersection& closestIntersection ) {

    // Don't check for intersection for large t
    float bound = std::numeric_limits<float>::max();

    closestIntersection.distance = bound;

    for (int i = 0; i < triangles.size(); i++) {
      vec4 v0 = triangles[i].v0;
      vec4 v1 = triangles[i].v1;
      vec4 v2 = triangles[i].v2;

      vec3 e1 = vec3(v1.x - v0.x ,v1.y - v0.y,v1.z - v0.z);
      vec3 e2 = vec3(v2.x - v0.x ,v2.y - v0.y,v2.z - v0.z);
      vec3 b = vec3(start.x-v0.x, start.y-v0.y, start.z-v0.z);

      //std::cout << to_string(b[0]) + " " + to_string(b[1]) + " " + to_string(b[2]) << '\n';

      // bottleneck
      // glm::mat3 A;
      vec3 dir3 = vec3(dir[0], dir[1], dir[2]);
      glm::mat3 A( -dir3, e1, e2 );

      vec3 x = glm::inverse( A ) * b;

      vec4 x4 = vec4(x[0], x[1], x[2], 1.0);

      bool check1 = x[1] >= 0;
      bool check2 = x[2] >= 0;
      bool check3 = (x[1] + x[2]) <= 1;
      bool check4 = x[0] >= 0;
      bool check5 = x[0] < bound;
      bool check6 = x[0] < closestIntersection.distance;

      if ( check1 && check2 && check3 && check4 && check5 && check6 ) {
        closestIntersection.position = x4;
        closestIntersection.distance = x[0];
        closestIntersection.triangleIndex = i;
      }
    }

    if ( closestIntersection.distance < bound ) {
        return true;
    }

    else
    {
      return false;
    }
  }
