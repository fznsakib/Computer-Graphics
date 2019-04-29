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

// Change to 1280 x 720 later
#define SCREEN_WIDTH 320
#define SCREEN_HEIGHT 256
#define FULLSCREEN_MODE false

// Done
// - Anti aliasing
// - Cramer's rule
// - Load sphere

// To do
// - Photon mapping
// - Convert point light to area
// - Fix rotation
// - Bump maps?

/* ----------------------------------------------------------------------------*/
/* STRUCTS                                                                   */
/* ----------------------------------------------------------------------------*/

struct Intersection {
   vec4 position;
   float distance;
   vec4 normal;
   vec3 colour;
   vec3 material;
   int triangleIndex;
   int sphereIndex;
};

struct Light {
  vec4 position;
  vec3 colour;
  float power;
};

struct Photon {
  vec4 position;
  vec4 direction;
  vec3 power;
  float phi;
  float theta;
  bool flag;   //flag used in kdtree
};

// phi = 255 * (atan2(dy,dx)+PI) / (2*PI)
// theta = 255 * acos(dx) / PI

/* ----------------------------------------------------------------------------*/
/* VARIABLES                                                                   */
/* ----------------------------------------------------------------------------*/

float focalLength = 256;
vec4 cameraPos( 0.0, 0.0, -3.0, 1.0);
vector<Light> lights;
float yaw = 0.0;
mat4 R(1.0f);

vector<Triangle> triangles;
vector<Sphere> spheres;

vector<Photon> causticPhotonMap;
vector<Photon> globalPhotonMap;
vector<Photon> volumePhotonMap;


/* ----------------------------------------------------------------------------*/
/* FUNCTIONS                                                                   */
/* ----------------------------------------------------------------------------*/


bool Update();
void Draw(screen* screen);
bool ClosestIntersection(vec4 start, vec4 dir, Intersection& closestIntersection );
vec3 DirectLight( const Intersection &i, Light light );
void PhotonEmission( const int noOfPhotons );
void TracePhoton( Photon& photon );
void GetReflectedDirection(const vec4& incident, const vec4& normal, vec4& reflected);
float RandomFloat(float min, float max);


int main( int argc, char* argv[] )
{

  screen *screen = InitializeSDL( SCREEN_WIDTH, SCREEN_HEIGHT, FULLSCREEN_MODE );

  // Initialise lights
  Light light;
  light.position = vec4( 0, -0.5, -0.7, 1.0 );
  light.colour = vec3( 14.f * vec3( 1, 1, 1 ));
  light.power = 10.0f;
  lights.push_back(light);

  // Initialise surfaces
  LoadTestModel( triangles, spheres );

  // Initialise photon map before scene rendering
  int noOfPhotons= 1000;
  PhotonEmission(noOfPhotons);

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
  // vector<Triangle> triangles;
  // vector<Sphere> spheres;
  //
  // LoadTestModel( triangles, spheres );


  // u and v are coordinates on the 2D screen
  for (int v = 0; v < SCREEN_HEIGHT; v++) {
    for (int u = 0; u < SCREEN_WIDTH; u++) {

      vec4 dir = (vec4((u - SCREEN_WIDTH/2), v - SCREEN_HEIGHT/2, focalLength, 1.0));
      /* Calculate the new direction vector after yaw rotation */
      dir = R * dir;

      pixelColour = vec3 (0.0f, 0.0f, 0.0f);
      bool validRay = false;

      // Trace multiple rays around each ray and average color for anti aliasing
      // Change to 4 x 4 later
      for (int i = -1; i <= 1; ++i) {
        for (int j = -1; j <= 1; ++j ) {
          // Multiply i and j by distance to surrounding rays being used
          float multiplier = 0.5;
          vec4 newDir = vec4(dir.x + (multiplier * i), dir.y + (multiplier * j), focalLength, 1.0);

          Intersection intersection;
          bool closestIntersection = ClosestIntersection(cameraPos, newDir, intersection);

          // If light intersects with object then draw it
          if (closestIntersection == true) {
            validRay = true;
            vec3 objectColor;

            // Check if intersection with triangle or sphere to get correct colour
            if (intersection.triangleIndex != -1) objectColor = triangles[intersection.triangleIndex].color;
            else objectColor = spheres[intersection.sphereIndex].color;

            for (int i = 0; i < lights.size(); i++) {
              pixelColour += DirectLight( intersection, lights[i] );
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

      float t;
      if (triangles[i].intersect(start3, dir3, t, closestIntersection.distance)) {
        // Position of intersection
        vec4 position = start + vec4(t * dir3, 0);

        // Check if distance to triangle is closest
        if (t < closestIntersection.distance) {
          closestIntersection.position = position;
          closestIntersection.distance = t;
          closestIntersection.normal = triangles[i].normal;
          closestIntersection.colour = triangles[i].color;
          closestIntersection.material = triangles[i].material;
          closestIntersection.triangleIndex = i;
          closestIntersection.sphereIndex = -1;
        }
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
          spheres[i].getNormal(vec3(position), closestIntersection.normal);
          closestIntersection.colour = spheres[i].color;
          closestIntersection.material = spheres[i].material;
          closestIntersection.triangleIndex = -1;
          closestIntersection.sphereIndex = i;
        }
      }
    }

    // Returns with intersections closest to origin of ray
    if ( closestIntersection.distance < bound ) {
        return true;
    }
    else {
      return false;
    }
  }


vec3 DirectLight( const Intersection &i, Light light ) {

  vec3 power, objectColor;
  vec4 normal;
  vec4 r = (light.position - i.position);
  float r_magnitude = sqrt(pow(r[0],2) + pow(r[1],2) + pow(r[2],2));

  vec4 direction = light.position - i.position;

  // Get normal and colour
  objectColor = i.colour;
  normal = i.normal;

  // Check for surface between triangle and light
  Intersection intersection;

  // Check if distance less than to light. Emit ray from a small distance off the surface
  // so that it doesn't intersect with itself
  if (ClosestIntersection(i.position + (normal * 0.00001f), direction, intersection)) {
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

// TODO
// - Trace photons

void PhotonEmission( const int noOfPhotons ) {
  int photonsEmitted = 0;
  float x, y, z;

  // Keep emitting photons until max number reached
  while (photonsEmitted < noOfPhotons) {
    Photon photon;

    // Rejection sampling to find random photon direction
    do {
      x = RandomFloat(-1.0, 1.0);
      y = RandomFloat(-1.0, 1.0);
      z = RandomFloat(-1.0, 1.0);
    } while ((x*x) + (y*y) + (z*z) > 1);

    // Change to square light later
    photon.position = vec4(lights[0].position.x, lights[0].position.y, lights[0].position.z, 1.0f);
    photon.direction = vec4(x, y, z, 1.0);
    photon.power = lights[0].colour / (float)noOfPhotons;

    // Trace photon from light position in photon direction dir
    TracePhoton(photon);

    photonsEmitted += 1;
  }
}


void TracePhoton(Photon &photon) {
  // Keep bouncing photon until it reaches diffuse object
  bool photonAbsorbed = false;
  Intersection intersection;

  while (photonAbsorbed == false) {
    // Call ClosestIntersection to find next object to intersect with
    ClosestIntersection(photon.position, photon.direction, intersection);

    // Use Russian Roulette on material properties to decide to photon action
    float random = RandomFloat(0.0f, 1.0f);

    // Check specularity for reflection occurring. Add photon for each intersection
    if (random < intersection.material[2]) {
      // Reflect photon with new position in new direction
      photon.position = intersection.position;

      globalPhotonMap.push_back(photon);

      GetReflectedDirection(photon.direction, intersection.normal, photon.direction);
    }

    // TODO - ADD FOR TRANSMISSION/REFRACTANCE
    // Check specularity + diffuseness for transmission occurring
    // else if (random < (intersection.material[1] + intersection.material[2])) {
    //   // Refract new photon
    //   photon.position = intersection.position;
    // }

    // Else absorb photon
    else {
      photon.position = intersection.position;

      globalPhotonMap.push_back(photon);

      photonAbsorbed = true;
    }

    // Decide which photon map to add photon to

  }

}


void GetReflectedDirection(const vec4& incident, const vec4& normal, vec4& reflected) {
	reflected = incident - normal * (2 * glm::dot(incident, normal));
}


float RandomFloat(float min, float max) {
  float num = ((float) rand()) / (float) RAND_MAX;

  float range = max - min;
  return (num*range) + min;
}
