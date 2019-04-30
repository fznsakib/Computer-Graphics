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

// Resolution - Change to 1280 x 720 later
#define SCREEN_WIDTH 320
#define SCREEN_HEIGHT 256
#define FULLSCREEN_MODE false

// Maximum reflection/refraction depth
#define MAX_RAY_DEPTH 5

// Done
// - Anti aliasing
// - Cramer's rule
// - Load sphere

// To do
// - Reflectance and Refraction
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
vec3 indirectLight = 0.5f * vec3(1, 1, 1);
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
bool ClosestIntersection(vec4 start, vec4 dir, Intersection& closestIntersection, int& depth );
vec3 DirectLight( const Intersection &i, Light light );
vec3 IndirectLight( const vec3 pixelColour, const vec3 objectColour );
void PhotonEmission( const int noOfPhotons );
void TracePhoton( Photon& photon );
void GetReflectedDirection(const vec4& incident, const vec4& normal, vec4& reflected);
void GetRefractedDirection(const vec4& incident, const vec4& normal, const float ior, vec4& refracted);
float RandomFloat(float min, float max);
float mix(const float &a, const float &b, const float &mix);


int main( int argc, char* argv[] )
{

  screen *screen = InitializeSDL( SCREEN_WIDTH, SCREEN_HEIGHT, FULLSCREEN_MODE );

  // Initialise lights
  Light light;
  light.position = vec4( 0, -0.5, -0.7, 1.0 );
  light.colour = vec3( 10.0f * vec3( 1, 1, 1 ));
  lights.push_back(light);

  // Initialise surfaces
  LoadTestModel( triangles, spheres );

  // Initialise photon map before scene rendering
  int noOfPhotons = 200;
  // PhotonEmission(noOfPhotons);

  std::cout << globalPhotonMap.size() << '\n';

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
          int depth = 0;
          bool closestIntersection = ClosestIntersection(cameraPos, newDir, intersection, depth);

          // If light intersects with object then draw it
          if (closestIntersection == true) {
            validRay = true;
            vec3 objectColour;

            // Check if intersection with triangle or sphere to get correct colour
            if (intersection.triangleIndex != -1) objectColour = triangles[intersection.triangleIndex].color;
            else objectColour = spheres[intersection.sphereIndex].color;

            // Direct illumination
            for (int i = 0; i < lights.size(); i++) {
              pixelColour += DirectLight( intersection, lights[i] );
            }

            // Indirect illumination for diffuse objects
            if (intersection.material[1] == 0.0f) {
              pixelColour = IndirectLight(pixelColour, objectColour);
            }
          }
        }
      }
      if (validRay == true) {
        // Average light for the multiple rays for anti aliasing
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
                          Intersection& closestIntersection,
                          int& depth ) {

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
          closestIntersection.sphereIndex = i;
          closestIntersection.triangleIndex = -1;
        }
      }
    }

    if (closestIntersection.distance == bound) {
      closestIntersection.colour = vec3(0.0f, 0.0f, 0.0f);
      return false;
    }

    // Reflection
    if (closestIntersection.material[1] > 0.0f) {
      Intersection reflection;

      // Calculate new reflected direction of ray
      vec4 reflectionDir = glm::reflect(dir, closestIntersection.normal);

      // Trace reflected ray
      if (ClosestIntersection(closestIntersection.position + (closestIntersection.normal * 0.00001f), reflectionDir, reflection, depth)) {
        vec3 colour = IndirectLight(closestIntersection.colour, reflection.colour);
        reflection.colour = colour;
      }

      closestIntersection.colour = reflection.colour;

      if (closestIntersection.colour == vec3(0.0f, 0.0f, 0.0f)) {
        return false;
      }
    }

    // REFLECTION AND REFRACTION
    // // If the normal and the view direction are not opposite to each other
    // // reverse the normal direction. That also means we are inside the sphere so set
    // // the inside bool to true. Finally reverse the sign of IdotN which we want
    // // positive.
    // bool inside = false;
    // if (glm::dot(dir, closestIntersection.normal) > 0) {
    //   closestIntersection.normal = -closestIntersection.normal;
    //   inside = true;
    // }
    //
    // float transmission = 1.0f - (closestIntersection.material[1] + closestIntersection.material[2]);
    // Intersection reflection;
    // Intersection refraction;
    //
    // // Check if reflective or refractive
    // if ((closestIntersection.material[2] > 0.0f || transmission > 0.0f) && depth < MAX_RAY_DEPTH) {
    //   float facingratio = -glm::dot(dir, closestIntersection.normal);
    //
    //   // Change the reflection/refraction mix value to tweak the effect
    //   float fresnelEffect = mix(pow(1 - facingratio, 3), 1, 0.1);
    //
    //   // Calculate new reflected direction of ray
    //   vec4 reflectionDir;
    //   GetReflectedDirection(dir, closestIntersection.normal, reflectionDir);
    //
    //   glm::normalize(reflectionDir);
    //
    //   // Trace reflected ray
    //   depth += 1;
    //   ClosestIntersection(closestIntersection.position + (closestIntersection.normal * 0.00001f), reflectionDir, reflection, depth);
    //
    //   // Check if refractive
    //   if (transmission > 0.0f) {
    //     // Compute index of refraction depending on whether light is inside or outside object
    //     float ior = 1.5;
    //     float eta = (inside) ? ior : 1 / ior;
    //
    //     // Calculate new refracted direction of ray
    //     vec4 refractionDir;
    //     GetRefractedDirection(dir, closestIntersection.normal, eta, refractionDir);
    //
    //     // Trace refracted ray
    //     depth += 1;
    //     ClosestIntersection(closestIntersection.position + (closestIntersection.normal * 0.00001f), refractionDir, refraction, depth);
    //   }
    //
    //   // Now compute the mix of reflection and refraction colours
    //   vec3 computedColour = ((reflection.colour * fresnelEffect) +
    //                         (refraction.colour * (1 - fresnelEffect) * transmission))
    //                         * closestIntersection.colour;
    //
    //   closestIntersection.colour += computedColour;
    // }

    return true;
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
  int depth = 0;

  // Check if distance less than to light. Emit ray from a small distance off the surface
  // so that it doesn't intersect with itself
  if (ClosestIntersection(i.position + (normal * 0.00001f), direction, intersection, depth)) {
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


vec3 IndirectLight( const vec3 pixelColour, const vec3 objectColour ) {
  return (pixelColour + (objectColour * indirectLight));
}


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
    int depth = 0;
    ClosestIntersection(photon.position, photon.direction, intersection, depth);

    // Use Russian Roulette on material properties to decide to photon action
    float random = RandomFloat(0.0f, 1.0f);

    // Add photon for each intersection
    // DIFFUSE REFLECTION
    if (random < intersection.material[0]) {
      // Reflect photon with new position in new direction
      photon.position = intersection.position;
      photon.power = photon.power * intersection.colour;

      globalPhotonMap.push_back(photon);

      // Need to change to random
      GetReflectedDirection(photon.direction, intersection.normal, photon.direction);
    }
    // SPECULAR REFLECTION
    else if (random < (intersection.material[1])) {
      // Reflect photon with new position in new direction
      photon.position = intersection.position;
      GetReflectedDirection(photon.direction, intersection.normal, photon.direction);

      // Don't store photon as this will be done by the raytracer
    }
    // ABSORPTION
    else {
      photon.position = intersection.position;
      photonAbsorbed = true;
    }
  }
}


void GetReflectedDirection(const vec4& incident, const vec4& normal, vec4& reflected) {
	reflected = incident - normal * (2 * glm::dot(incident, normal));
}


void GetRefractedDirection(const vec4& incident, const vec4& normal, const float ior, vec4& refracted) {
  float cosi = -(glm::dot(normal, incident));
  float k = (1 - ior) * ior * (1 - cosi * cosi);

  refracted = (incident * ior) + normal * (ior *  cosi - sqrt(k));
  glm::normalize(refracted);
}


float RandomFloat(float min, float max) {
  float num = ((float) rand()) / (float) RAND_MAX;

  float range = max - min;
  return (num*range) + min;
}


float mix(const float &a, const float &b, const float &mix)
{
    return b * mix + a * (1 - mix);
}
