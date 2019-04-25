#include <iostream>
#include <glm/glm.hpp>
#include <SDL.h>
#include <stdint.h>
#include <limits.h>
#include <math.h>

#include "SDLauxiliary.h"
#include "TestModelH.h"

using namespace std;
using glm::vec3;
using glm::mat3;
using glm::vec4;
using glm::mat4;

SDL_Event event;

#define SCREEN_WIDTH 320
#define SCREEN_HEIGHT 256
#define FULLSCREEN_MODE false

// Done
// - Anti aliasing
// - Cramer's rule
// - Load sphere

// 

// To do
// - Clean up code
// - Photon mapping
// - Fix rotation
// - Bump maps?

/* ----------------------------------------------------------------------------*/
/* STRUCTS                                                                   */
/* ----------------------------------------------------------------------------*/

struct Intersection {
       vec4 position;
       float distance;
       int triangleIndex;
       int sphereIndex;
};

struct Light {
      vec4 position;
      vec3 colour;
};

/* ----------------------------------------------------------------------------*/
/* VARIABLES                                                                   */
/* ----------------------------------------------------------------------------*/

float focalLength = 256;
vec4 cameraPos( 0.0, 0.0, -3.0, 1.0);
vector<Light> lights;
float yaw = 0.0;
mat4 R(1.0f);


/* ----------------------------------------------------------------------------*/
/* FUNCTIONS                                                                   */
/* ----------------------------------------------------------------------------*/


bool Update();
void Draw(screen* screen);
bool ClosestIntersection(vec4 start, vec4 dir,
                         const vector<Triangle>& triangles,
                         const vector<Sphere>& spheres,
                         Intersection& closestIntersection );
vec3 DirectLight( const Intersection &i,
                  const vector<Triangle>& triangles,
                  const vector<Sphere>& spheres,
                  Light light );


int main( int argc, char* argv[] )
{

  screen *screen = InitializeSDL( SCREEN_WIDTH, SCREEN_HEIGHT, FULLSCREEN_MODE );

  // Initialise lights
  Light light1;
  light1.position = vec4( 0, -0.5, -0.7, 1.0 );
  light1.colour = vec3( 14.f * vec3( 1, 1, 1 ));
  lights.push_back(light1);

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

  vec3 black(0.0,0.0,0.0);
  vec3 pixelColour = vec3( 1.0, 1.0, 1.0 );
  vec3 indirectLight = 0.5f * vec3(1, 1, 1);

  // Initialise triangles and spheres
  vector<Triangle> triangles;
  vector<Sphere> spheres;

  LoadTestModel( triangles, spheres );

  // for (int i = 0; i < triangles.size(); i++) {
  //   objects.push_back(std::unique_ptr<Object> (new Triangle(triangles[i].v0, triangles[i].v1, triangles[i].v2, triangles[i].color)));
  // }

  // u and v are coordinates on the 2D screen
  for (int v = 0; v < SCREEN_HEIGHT; v++) {
    for (int u = 0; u < SCREEN_WIDTH; u++) {

      vec4 dir = (vec4((u - SCREEN_WIDTH/2), v - SCREEN_HEIGHT/2, focalLength, 1.0));
      /* Calculate the new direction vector after yaw rotation */
      dir = R * dir;

      pixelColour = vec3 (0.0f, 0.0f, 0.0f);
      bool validRay = false;

      // Trace multiple rays around each ray and average color for anti aliasing
      for (int i = -1; i <= 1; ++i) {
        for (int j = -1; j <= 1; ++j ) {
          float multiplier = 0.5;
          vec4 newDir = vec4(dir.x + (multiplier * i), dir.y + (multiplier * j), focalLength, 1.0);

          Intersection intersection;
          bool closestIntersection = ClosestIntersection(cameraPos, newDir, triangles, spheres, intersection);

          // If light intersects with object then draw it
          if (closestIntersection == true) {
            validRay = true;
            vec3 objectColor;

            // Check if intersection with triangle or sphere to get correct colour
            if (intersection.triangleIndex != -1) objectColor = triangles[intersection.triangleIndex].color;
            else objectColor = spheres[intersection.sphereIndex].color;

            for (int i = 0; i < lights.size(); i++) {
              pixelColour += DirectLight( intersection, triangles, spheres, lights[i] );
            }

            // Take into account indirect illumination
            pixelColour = pixelColour + (objectColor * indirectLight);
          }
        }
      }
      if (validRay == true) {
        vec3 averageLight = pixelColour/9.0f;
        PutPixelSDL(screen, u, v, averageLight);
      }
      else {
        PutPixelSDL(screen, u, v, black);
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
          // LIGHT MOVEMENT
          case SDLK_w:
            lights[0].position += vec4(0, 0, 0.1, 0);
          break;
          case SDLK_s:
            lights[0].position += vec4(0, 0, -0.1, 0);
          break;
          case SDLK_a:
            lights[0].position += vec4(-0.1, 0, 0, 0);
          break;
          case SDLK_d:
            lights[0].position += vec4(0.1, 0, 0, 0);
            // std::cout << lights[0].position[0] << '\n';
          break;
          case SDLK_q:
            lights[0].position += vec4(0, -0.1, 0, 0);
          break;
          case SDLK_e:
            lights[0].position += vec4(0, 0.1, 0, 0);
          break;

          // CAMERA MOVEMENT
	        case SDLK_UP:
		        cameraPos += vec4(0, 0, 0.1, 0);
            cout << "KEYPRESS UP." <<endl;
		      break;
	        case SDLK_DOWN:
		        cameraPos += vec4(0, 0, -0.1, 0);
            cout << "KEYPRESS DOWN." <<endl;
		      break;
	        case SDLK_LEFT:
		        cameraPos += vec4(-0.1, 0, 0, 0);
            cout << "KEYPRESS <-." <<endl;
		      break;
	        case SDLK_RIGHT:
		        cameraPos += vec4(0.1, 0, 0, 0);
            cout << "KEYPRESS ->." <<endl;
		      break;
          // CAMERA ROTATION
          case SDLK_n:
            cout << "KEYPRESS a." <<endl;
            yaw -= 0.174533;
            R[0][0] = cos(yaw);   R[0][1] = 0;   R[0][2] = -sin(yaw);
            R[1][0] = 0;          R[1][1] = 1;   R[1][2] = 0;
            R[2][0] = sin(yaw);   R[2][1] = 0;   R[2][2] = cos(yaw);
          break;
          case SDLK_m:
            yaw += 0.174533;
            R[0][0] = cos(yaw);   R[0][1] = 0;   R[0][2] = -sin(yaw);
            R[1][0] = 0;          R[1][1] = 1;   R[1][2] = 0;
            R[2][0] = sin(yaw);   R[2][1] = 0;   R[2][2] = cos(yaw);
          break;
          // FOCAL LENGTH CHANGE
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


bool ClosestIntersection( vec4 start, vec4 dir,
                          const vector<Triangle>& triangles,
                          const vector<Sphere>& spheres,
                          Intersection& closestIntersection ) {

    // Don't check for intersection for large t
    float bound = std::numeric_limits<float>::max();

    closestIntersection.distance = bound;
    vec3 start3 = vec3(start[0], start[1], start[2]);
    vec3 dir3 = vec3(dir[0], dir[1], dir[2]);

    ///////////////
    // TRIANGLES //
    ///////////////
    for (int i = 0; i < triangles.size(); i++) {
      vec4 v0 = triangles[i].v0;
      vec4 v1 = triangles[i].v1;
      vec4 v2 = triangles[i].v2;

      vec3 e1 = vec3(v1.x - v0.x ,v1.y - v0.y,v1.z - v0.z);
      vec3 e2 = vec3(v2.x - v0.x ,v2.y - v0.y,v2.z - v0.z);
      // vec3 b = vec3(start.x-v0.x, start.y-v0.y, start.z-v0.z);

      // Currently: 55000ms without -O3, 600ms with -O3
      // With Cramer's rule: 47000ms without -O3, 430ms with -O3
      glm::mat3 system( -dir3, e1, e2 );

      // Use Cramer's rule to produce just t (x[0]). If t is negative, then
      // return false. Else, continue with calculation

      // Using the linear equation:
      // (-dir, e1, e2) (t, u, v) = s - v0
      vec4 solutions = start - v0;
      vec3 solutions3 = vec3(solutions[0], solutions[1], solutions[2]);

      // Solving t in using Cramer's rule:
      // -dir[0]t + e1[0]u + e2[0]v = solutions[0]
      // -dir[1]t + e1[1]u + e2[1]v = solutions[0]
      // -dir[2]t + e1[2]u + e2[2]v = solutions[0]

      // Where t = det(system_t)/det(system)
      glm::mat3 system_t( solutions3, e1, e2 );
      float t = glm::determinant(system_t)/glm::determinant(system);
      float distance = t * glm::length(dir3);


      // If the distance is negative, the intersection does not occur, so carry on to the next one.
      if (distance < 0.0f) continue;
      // Do distance checks earlier to move on quicker if possible
      else if (distance >= closestIntersection.distance || distance > bound) continue;


      // Similar to t, calculate u and v
      glm::mat3 system_u( -dir3, solutions3, e2 );
      float u = glm::determinant(system_u)/glm::determinant(system);

      glm::mat3 system_v( -dir3, e1, solutions3 );
      float v = glm::determinant(system_v)/glm::determinant(system);

      // x = (t, u, v)
      // vec3 x = glm::inverse( system ) * b;

      vec4 position = start + vec4(t * dir3, 0);

      bool check = (u >= 0) && (v >= 0) && ((u + v) <= 1);

      if (check) {
        closestIntersection.position = position;
        closestIntersection.distance = distance;
        closestIntersection.triangleIndex = i;
        closestIntersection.sphereIndex = -1;
      }
    }

    /////////////
    // SPHERES //
    /////////////
    for (int i = 0; i < spheres.size(); i++) {
      float t;
      if (spheres[i].intersect(start3, dir3, t)) {
        // Position of intersection
        vec4 position = start + vec4(t * dir3, 0);

        // Check if distance to sphere is closest
        if (t < closestIntersection.distance) {
          closestIntersection.position = position;
          closestIntersection.distance = t; //Might be t * dir3 as above
          closestIntersection.triangleIndex = -1;
          closestIntersection.sphereIndex = i;
        }
      }
    }

    if ( closestIntersection.distance < bound ) {
        return true;
    }
    else {
      return false;
    }
  }


vec3 DirectLight( const Intersection &i, const vector<Triangle>& triangles, const vector<Sphere>& spheres, Light light ) {

  vec3 power, objectColor;
  vec4 normal;
  vec4 r = (light.position - i.position);
  float r_magnitude = sqrt(pow(r[0],2) + pow(r[1],2) + pow(r[2],2));

  vec4 direction = light.position - i.position;

  // Check if object is triangle of sphere for normakl and colour
  // Get normal and colour
  if (i.triangleIndex != -1) {
    objectColor = triangles[i.triangleIndex].color;
    normal = triangles[i.triangleIndex].normal;
  }
  else {
    objectColor = spheres[i.sphereIndex].color;
    vec3 normal3;
    spheres[i.sphereIndex].getNormal(vec3(i.position), normal3);

    normal = vec4(normal3.x, normal3.y, normal3.z, 0);
  }

  // Check for surface between triangle and light
  Intersection intersection;

  // Check if distance less than to light. Emit ray from a small distance off the surface
  // so that it doesn't intersect with itself
  if (ClosestIntersection(i.position + (normal * 0.00001f), direction, triangles, spheres, intersection)) {
    if (intersection.distance < r_magnitude) {
      return vec3(0.0,0.0,0.0);
    }
  }

  vec3 normalisedDirection = glm::normalize(vec3(direction));
  vec3 normal3 = vec3(normal);

  float a = glm::dot(normalisedDirection, normal3);
  float b = 4 * M_PI;

  float surfaceArea = (b * pow(r_magnitude,2));

  // If light does not hit triangle
  if (a <= 0) a = 0.f;

  // Return power of direct illumination
  power = (objectColor * light.colour * a)/surfaceArea;

  return power;
}
