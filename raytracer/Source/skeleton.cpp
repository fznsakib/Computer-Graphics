#include <iostream>
#include <glm/glm.hpp>
#include <SDL.h>
#include <stdint.h>
#include <limits.h>
#include <math.h>
#include <algorithm>
#include <functional>
#include <array>
#include <iostream>

#include "SDLauxiliary.h"
#include "TestModelH.h"

using namespace std;
using glm::vec3;
using glm::mat3;
using glm::vec4;
using glm::mat4;

SDL_Event event;

// Resolution
// Default - 320 x 256 / focalLength = 256 / z = -3.0
// Medium  - 800 x 640
// Change to 1280 x 720 / focalLength = 1024 / z = -4.2
#define SCREEN_WIDTH 320
#define SCREEN_HEIGHT 256
#define FULLSCREEN_MODE false


// Maximum reflection/refraction depth
#define MAX_RAY_DEPTH 20

// Done
// - Anti aliasing
// - Cramer's rule
// - Load sphere
// - Convert point light to area
// - Reflection
// - Refraction
// - Photon Mapping

// To do
// - Photon Mapping with reflection and refraction
// - Colour bleeding in photon mapping
// - Build kd tree??
// - Procedural floors/walls

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
  vec3 position;
  vec3 direction;
  vec3 power;
  float distance;
  bool flag;   //flag used in kdtree
};


struct KDTree {
  Photon* node;
  KDTree* left;
  KDTree* right;
};


/* ----------------------------------------------------------------------------*/
/* VARIABLES                                                                   */
/* ----------------------------------------------------------------------------*/

const float TRANSLATION = 0.25f;
const float ROTATION = 0.125f;

float focalLength = 256;
vec4 cameraPos( 0.0, 0.0, -3.0, 1.0);
vector<Light> lights;
float yaw = 0.0;
mat4 R(1.0f);
vec3 black(0.0f,0.0f,0.0f);
vec3 indirectLight = 0.5f * vec3(1, 1, 1);


vector<Triangle> triangles;
vector<Sphere> spheres;

vector<Photon> globalPhotonMap;

vector<Photon*> globalPhotonPointers;

KDTree* globalPhotonTree;

int currentMaxDimension;
int radianceCount = 500;
int noNode = -1;
float searchRadius = 0.1f;
int noOfPhotons = 20000;

KDTree* endTree = (KDTree*)malloc(1 * sizeof(KDTree));
Photon* endPhoton = (Photon*)malloc(1 * sizeof(Photon));


/* ----------------------------------------------------------------------------*/
/* FUNCTIONS                                                                   */
/* ----------------------------------------------------------------------------*/


bool Update();
void Draw(screen* screen);
vec3 CastRay( const vec4 &position, const vec4 &dir, const int& depth );
bool ClosestIntersection(vec4 start, vec4 dir, Intersection& closestIntersection );
vec3 DirectLight( const Intersection &i, Light light );
vec3 IndirectLight( const vec3 pixelColour, const vec3 objectColour, const vec3 light );
void PhotonEmission( const int noOfPhotons );
void TracePhoton( Photon& photon );
void GetReflectedDirection(const vec4& incident, const vec4& normal, vec4& reflected);
void GetRefractedDirection(const vec4& incident, const vec4& normal, const float ior, vec4& refracted);
void GetFresnel(const vec4& incident, const vec4& normal, const float& ior, float& kr);
vec4 GenerateRandomDirection();
KDTree* BalanceTree( vector<Photon*> photonPointers );
void LocatePhotons( Photon* rootNode, vector<Photon*> maxHeap, Intersection intersection, float& maxSearchDistance );
int MaxDimension( vector<Photon> photons );
int GetCurrentMaxDimension();
float RandomFloat(float min, float max);
float Mix(const float &a, const float &b, const float &mix);
void MaxHeapSort( vector<Photon*>& maxHeap, const Intersection intersection);
void Heapify(vector<Photon*>& maxHeap, int heapSize, int i, vec3 position, vec3 normal);
vec3 RadianceEstimate( vec4 intersectionPosition );


int main( int argc, char* argv[] ) {

  screen *screen = InitializeSDL( SCREEN_WIDTH, SCREEN_HEIGHT, FULLSCREEN_MODE );

  // Initialise lights
  Light light;
  light.position = vec4( 0, -0.2, -0.7, 1.0 );
  light.colour = vec3( 500.0f * vec3( 1, 1, 1 ));
  lights.push_back(light);

  // Initialise surfaces
  LoadTestModel( triangles, spheres );

  // Initialise photon map before scene rendering
  // int noOfPhotons = 500000;
  PhotonEmission(noOfPhotons);

  // Populate vector with pointers to photons
  for (int i = 0; i < globalPhotonMap.size(); i++) {
    globalPhotonPointers.push_back(&globalPhotonMap[i]);
  }

  // Create and balance photon kd tree
  globalPhotonTree = BalanceTree(globalPhotonPointers);


  while ( Update())
    {
      Draw(screen);
      SDL_Renderframe(screen);
    }

  SDL_SaveImage( screen, "screenshot.bmp" );

  KillSDL(screen);
  return 0;
}

vec3 RadianceEstimate( vec4 intersectionPosition ) {
  // Locate nearest photons
  vector<Photon> nearestPhotons;
  vec3 totalRadiance;

  for (int i = 0; i < globalPhotonMap.size(); i++) {
    float distance = glm::distance(vec3(intersectionPosition), globalPhotonMap[i].position);
    if (distance < searchRadius) {
      globalPhotonMap[i].distance = distance;
      nearestPhotons.push_back(globalPhotonMap[i]);
    }
  }

  // // Sort photons by distance, take closest n photons
  sort(nearestPhotons.begin(), nearestPhotons.end(), [](const Photon lhs, const Photon rhs) {
    return lhs.distance < rhs.distance;
  });

  // Add up radiance from closest photons
  if (nearestPhotons.size() == 0) {
    totalRadiance = vec3(0.0f, 0.0f, 0.0f);
  }
  else if (nearestPhotons.size() > radianceCount || nearestPhotons.size() == radianceCount) {
    for (int i = 0; i < radianceCount; i++){
      totalRadiance += nearestPhotons[i].power;
    }
    // totalRadiance = totalRadiance / float((M_PI * (searchRadius * searchRadius)));
  }
  else {
    for (int i = 0; i < nearestPhotons.size(); i++) {
      totalRadiance += nearestPhotons[i].power;
    }
    // totalRadiance = totalRadiance / float((M_PI * (searchRadius * searchRadius)));
  }

  return totalRadiance;

  // vector<Photon*> nearestPhotons;
  // nearestPhotons.reserve(radianceCount);

  // LocatePhotons( globalPhotonTree, nearestPhotons, intersection, searchRadius );

  // Calculate radiance estimate
  // for (int i = 0; i < nearestPhotons.size(); i++ ) {
    // float distance = glm::dot(intersection.normal, (intersection.position - nearestPhotons[i]->position));
    // std::cout << "distance of photon no " << i << ": " << distance << '\n';
  //   totalRadiance += nearestPhotons[i]->power;
  // }
}

/*Place your drawing here*/
void Draw(screen* screen) {

  /* Clear buffer */
  memset(screen->buffer, 0, screen->height*screen->width*sizeof(uint32_t));

  // u and v are coordinates on the 2D screen
  for (int v = 0; v < SCREEN_HEIGHT; v++) {
    for (int u = 0; u < SCREEN_WIDTH; u++) {
      vec4 dir = (vec4((u - SCREEN_WIDTH/2), v - SCREEN_HEIGHT/2, focalLength, 1.0));
      // Calculate the new direction vector after yaw rotation
      dir = R * dir;
      glm::normalize(dir);


      // Apply global photon mapping
      Intersection intersection;
      bool lightCheck = false;
      vec3 totalRadiance = vec3(0.0f, 0.0f, 0.0f);

      if (ClosestIntersection(cameraPos, dir, intersection)) {
        // RENDER LIGHT EMITTER
        // lightCheck = intersection.position.x >  -0.2f  &&
        //              intersection.position.x <   0.2f  &&
        //              intersection.position.y >= -1.0f  &&
        //              intersection.position.y <  -0.95f &&
        //              intersection.position.z >  -0.8f  &&
        //              intersection.position.z <  -0.4f;
        //
        // if (lightCheck) {
        //   PutPixelSDL(screen, u, v, vec3(1.0f, 1.0f, 1.0f));
        //   continue;
        // }

        // RENDER SPECULAR OBJECTS
        if (intersection.material[1] > 0.0f || intersection.material[2] > 0.0f) {
          int depth = 0;
          totalRadiance = CastRay(cameraPos, dir, depth);
          PutPixelSDL(screen, u, v, totalRadiance);
          continue;
        }

        // RENDER GLOBAL PHOTON MAP
        totalRadiance = RadianceEstimate(intersection.position);

      }

      PutPixelSDL(screen, u, v, totalRadiance);
    }
  }
}


      // // OLD RENDERING METHOD
      // vec3 pixelColour = vec3(0.0f, 0.0f, 0.0f);
      // bool validRay = false;
      //
      // // Trace multiple rays around each ray and average color for anti aliasing
      // // Change to 4 x 4 later
      // for (int i = -2; i <= 2; ++i) {
      //   for (int j = -2; j <= 2; ++j ) {
      //     // Multiply i and j by distance to surrounding rays being used
      //     float multiplier = 0.25;
      //     vec4 newDir = vec4(dir.x + (multiplier * i), dir.y + (multiplier * j), focalLength, 1.0);
      //
      //     Intersection intersection;
      //     int depth = 0;
      //     pixelColour += CastRay(cameraPos, newDir, depth);
      //
      //     validRay = true;
      //
      //     // bool closestIntersection = ClosestIntersection(cameraPos, newDir, intersection);
      //     //
      //     // // If light intersects with object then draw it
      //     // if (closestIntersection == true) {
      //     //   validRay = true;
      //     //   vec3 objectColour = intersection.colour;
      //     //
      //     //   // Direct illumination
      //     //   for (int i = 0; i < lights.size(); i++) {
      //     //     if (DirectLight( intersection, lights[i] ) != vec3(0.0f, 0.0f, 0.0f)) {
      //     //       pixelColour += DirectLight( intersection, lights[i] );
      //     //     }
      //     //   }
      //     //
      //     //   // Indirect illumination for diffuse objects
      //     //   if (intersection.material[1] == 0.0f) {
      //     //     pixelColour = IndirectLight(pixelColour, objectColour, indirectLight);
      //     //   }
      //     // }
      //   }
      // }
      // if (validRay == true) {
        // Average light for the multiple rays for anti aliasing
        // vec3 averageLight = pixelColour/16.0f;
        // PutPixelSDL(screen, u, v, averageLight);
      // }
      // else {
      //   PutPixelSDL(screen, u, v, black);
      // }
//     }
//   }
// }


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
            lights[0].position += vec4(0, 0, TRANSLATION, 0);
          break;
          case SDLK_s:
            lights[0].position += vec4(0, 0, -TRANSLATION, 0);
          break;
          case SDLK_a:
            lights[0].position += vec4(-TRANSLATION, 0, 0, 0);
          break;
          case SDLK_d:
            lights[0].position += vec4(TRANSLATION, 0, 0, 0);
          break;
          case SDLK_q:
            lights[0].position += vec4(0, -TRANSLATION, 0, 0);
          break;
          case SDLK_e:
            lights[0].position += vec4(0, TRANSLATION, 0, 0);
          break;

          // CAMERA MOVEMENT
	        case SDLK_UP:
		        cameraPos += vec4(0, 0, TRANSLATION, 0);
		      break;
	        case SDLK_DOWN:
		        cameraPos += vec4(0, 0, -TRANSLATION, 0);
		      break;
	        case SDLK_LEFT:
		        cameraPos += vec4(-TRANSLATION, 0, 0, 0);
		      break;
	        case SDLK_RIGHT:
		        cameraPos += vec4(TRANSLATION, 0, 0, 0);
		      break;

          // CAMERA ROTATION
          case SDLK_n:
            cout << "KEYPRESS n" <<endl;
            yaw -= ROTATION;
            R[0][0] = cos(yaw);   R[0][1] = 0;   R[0][2] = sin(yaw);
            R[1][0] = 0;          R[1][1] = 1;   R[1][2] = 0;
            R[2][0] = -sin(yaw);   R[2][1] = 0;   R[2][2] = cos(yaw);
          break;
          case SDLK_m:
            cout << "KEYPRESS m" <<endl;
            yaw += ROTATION;
            R[0][0] = cos(yaw);   R[0][1] = 0;   R[0][2] = sin(yaw);
            R[1][0] = 0;          R[1][1] = 1;   R[1][2] = 0;
            R[2][0] = -sin(yaw);   R[2][1] = 0;   R[2][2] = cos(yaw);
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


vec3 CastRay( const vec4 &position, const vec4 &dir, const int& depth ) {
  if (depth > MAX_RAY_DEPTH) return black;

  vec3 hitColour = vec3(0.0f, 0.0f, 0.0f);
  Intersection intersection;

  // Check if ray intersects with anything
  if (ClosestIntersection(position, dir, intersection)) {
    // Decide what to do depending on material

    // DIFFUSE
    if(intersection.material[0] == 0.5f) {
      // vec3 illumination = vec3(0.0f, 0.0f, 0.0f);

      // Compute direct light
      // for (int i = 0; i < lights.size(); i++) {
      //   if (DirectLight( intersection, lights[i] ) != vec3(0.0f, 0.0f, 0.0f)) {
      //     illumination += DirectLight( intersection, lights[i] );
      //   }
      // }

      // Add direct light
      // hitColour += illumination;

      // Add indirect light
      // hitColour = IndirectLight(hitColour, intersection.colour, indirectLight);

      hitColour = RadianceEstimate(intersection.position);
    }
    // ONLY REFLECTION
    else if(intersection.material[1] > 0.0f && intersection.material[2] == 0.0f) {

      // Compute the new reflected direction
      vec4 reflectedDir;
      GetReflectedDirection(dir, intersection.normal, reflectedDir);
      glm::normalize(reflectedDir);

      hitColour += 0.8f * CastRay(intersection.position + (intersection.normal * 0.00001f), reflectedDir, depth + 1);
    }
    // REFLECTION AND REFRACTION
    else if(intersection.material[1] > 0.0f && intersection.material[2] > 0.0f) {
      vec3 reflectionColour = vec3(0.0f, 0.0f, 0.0f);
      vec3 refractionColour = vec3(0.0f, 0.0f, 0.0f);

      // Index of refraction for glass
      float ior = 1.2;

      // Compute fresnel
      // float kr = 1.0f;
      // GetFresnel(dir, intersection.normal, ior, kr);

      bool outside = (glm::dot(dir, intersection.normal) < 0);
      vec4 offset = 0.00001f * intersection.normal;

      // Compute refraction ray if it is not a case of total internal reflection
      // if (kr < 1) {
        vec4 refractedDir;
        GetRefractedDirection(dir, intersection.normal, ior, refractedDir);
        glm::normalize(refractedDir);
        vec4 refractionRayOrigin = outside ? intersection.position - offset : intersection.position + offset;
        // Recursively cast refracted ray until depth is reached
        refractionColour = CastRay(refractionRayOrigin, refractedDir, depth + 1);
      // }


      vec4 reflectedDir;
      GetReflectedDirection(dir, intersection.normal, reflectedDir);
      glm::normalize(reflectedDir);
      vec4 reflectionRayOrigin = outside ? intersection.position + offset : intersection.position - offset;
      reflectionColour = CastRay(reflectionRayOrigin, reflectedDir, depth + 1);

      // Mix the relection and refraction to get one colour
      // hitColour += (reflectionColour * kr) + (refractionColour * (1 - kr));
      // hitColour += refractionColour * (1 - kr);
      hitColour += refractionColour;

    }
  }
  // No intersection, so return black
  else {
    hitColour = black;
  }

  return hitColour;
}


bool ClosestIntersection( vec4 start, vec4 dir,
                          Intersection& closestIntersection) {

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

    // If no intersection then return black
    if (closestIntersection.distance == bound) {
      closestIntersection.colour = vec3(0.0f, 0.0f, 0.0f);
      return false;
    }
    else {
      return true;
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
      return vec3(0.0f,0.0f,0.0f);
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


vec3 IndirectLight( const vec3 pixelColour, const vec3 objectColour, const vec3 indirectLight ) {
  return (pixelColour + (objectColour * indirectLight));
}


void PhotonEmission( const int noOfPhotons ) {
  int photonsEmitted = 0;
  float x, y, z;

  // Keep emitting photons until max number reached
  while (photonsEmitted < noOfPhotons) {
    Photon photon;

    // Using square light
    photon.position = vec3(RandomFloat(-0.1f, 0.1f), -0.95f, RandomFloat(-0.2f, 0.0f));
    float u = RandomFloat(0.0f, 1.0f);
    float v = 2 * M_PI * RandomFloat(0.0f, 1.0f);
    photon.direction = vec3((cos(v) * sqrt(u)), sqrt(1- u), (sin(v) * sqrt(u)));


    photon.power = lights[0].colour / (float)noOfPhotons;

    float random = RandomFloat(0.0f, 1.0f);
    if (random < 0.5f) globalPhotonMap.push_back(photon);


    // Randomly choose RGB colour for colour bleeding photon map
    float randomColourPicker = rand() % 3;

    if (randomColourPicker == 0)      photon.power = vec3(1.0f, 0.0f, 0.0f);
    else if (randomColourPicker == 1) photon.power = vec3(0.0f, 1.0f, 0.0f);
    else                              photon.power = vec3(0.0f, 0.0f, 1.0f);


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

    if (!ClosestIntersection(vec4(photon.position, 1.0), vec4(photon.direction, 1.0), intersection)) {
      break;
    }

    // Use Russian Roulette on material properties to decide to photon action
    float random = RandomFloat(0.0f, 1.0f);

    // Add photon for each intersection
    // DIFFUSE REFLECTION
    if (random < intersection.material[0]) {
      // Reflect photon with new position in new direction
      photon.position = vec3(intersection.position);
      photon.power = photon.power * intersection.colour;

      // Absorb depending on colour
      globalPhotonMap.push_back(photon);

      // Reflect in random direction
      vec4 randomDirection = GenerateRandomDirection();
      Intersection intersectionTest;

      while (!ClosestIntersection(vec4(photon.position, 1.0f), randomDirection, intersectionTest)) {
        randomDirection = GenerateRandomDirection();
      }

      photon.direction = vec3(randomDirection);
    }
    // SPECULAR REFLECTION
    else if (random < (intersection.material[0] + intersection.material[1])) {
      // Reflect photon with new position in new direction
      photon.position = vec3(intersection.position);
      vec4 reflectedDir;
      vec4 incident = vec4(photon.direction.x, photon.direction.y, photon.direction.z, 1.0f);
      GetReflectedDirection(incident, intersection.normal, reflectedDir);
      photon.direction = vec3(reflectedDir);
      // Don't store photon as this will be done by the raytracer
    }
    // TRANSMISSION
    else if (random < (intersection.material[0] + intersection.material[1] + intersection.material[2])) {

      bool outside = (glm::dot(photon.direction, vec3(intersection.normal)) < 0);
      vec3 offset = 0.00001f * vec3(intersection.normal);

      vec4 refractedDir;
      vec4 incident = vec4(photon.direction.x, photon.direction.y, photon.direction.z, 1.0f);
      GetRefractedDirection(incident, intersection.normal, 1.5, refractedDir);

      vec3 refractionRayOrigin = outside ? vec3(intersection.position) - offset : vec3(intersection.position) + offset;

      photon.position = refractionRayOrigin;
      photon.direction = vec3(refractedDir);

    }
    // ABSORPTION
    else {
      photon.position = vec3(intersection.position);
      photonAbsorbed = true;

      // globalPhotonMap.push_back(photon);
    }
  }
}


void GetReflectedDirection(const vec4& incident, const vec4& normal, vec4& reflected) {
	reflected = incident - normal * (2 * glm::dot(incident, normal));
}


void GetRefractedDirection(const vec4& incident, const vec4& normal, const float ior, vec4& refracted) {

  float cosi = glm::clamp(-1.0f, 1.0f, glm::dot(normal, incident));
  // float cosi = -(glm::dot(normal, incident));

  // Refraction index for outside and inside the object
  float etai = 1, etat = ior;

  vec4 n = normal;

  if (cosi < 0) {
    cosi = -cosi;
  }
  else {
    std::swap(etai, etat);
    n = -normal;
  }

  float eta = etai / etat;
  float k = 1 - eta * eta * (1 - cosi * cosi);
  // return k < 0 ? 0 : eta * incident + (eta * cosi - sqrtf(k)) * n;

  eta = 2.0f - eta;
  cosi = -(glm::dot(normal, incident));

  // If angle greater than critical (k is negative) angle so there is not refraction
  if (k < 0) refracted = vec4(0.0f, 0.0f, 0.0f, 0.0f);
  // else       refracted = eta * incident + (eta * cosi - sqrtf(k)) * n;
  else refracted = (incident * eta - normal * (-cosi + eta * cosi));

  glm::normalize(refracted);
}


void GetFresnel(const vec4& incident, const vec4& normal, const float& ior, float& kr) {

  float cosi = glm::clamp(-1.0f, 1.0f, glm::dot(incident, normal));
  float etai = 1, etat = ior;

  if (cosi > 0) {
    std::swap(etai, etat);
  }

  // Compute sini using Snell's law
  float sint = etai / etat * sqrtf(std::max(0.f, 1 - cosi * cosi));
  // Total internal reflection
  if (sint >= 1) {
      kr = 1;
  }
  else {
      float cost = sqrtf(std::max(0.f, 1 - sint * sint));
      cosi = fabsf(cosi);
      float Rs = ((etat * cosi) - (etai * cost)) / ((etat * cosi) + (etai * cost));
      float Rp = ((etai * cosi) - (etat * cost)) / ((etai * cosi) + (etat * cost));
      kr = (Rs * Rs + Rp * Rp) / 2;
  }

  // Due to conservation of energy, amout of refraction is:
  // kt = 1 - kr;
}


vec4 GenerateRandomDirection() {
  float x, y, z;
  do {
    x = RandomFloat(-1.0, 1.0);
    y = RandomFloat(-1.0, 1.0);
    z = RandomFloat(-1.0, 1.0);
  } while ((x*x) + (y*y) + (z*z) > 1);

  return vec4(x, y, z, 1.0);
}


float RandomFloat(float min, float max) {
  float num = ((float) rand()) / (float) RAND_MAX;

  float range = max - min;
  return (num*range) + min;
}


float Mix(const float &a, const float &b, const float &mix) {
    return b * mix + a * (1 - mix);
}

int GetCurrentMaxDimension() {
  return currentMaxDimension;
}

int MaxDimension( vector<Photon*> photons ) {
  float minX =  2.0f, minY =  2.0f, minZ =  2.0f;
  float maxX = -2.0f, maxY = -2.0f, maxZ = -2.0f;
  // Find the dimension where distance is max
  for (int i = 0; i < photons.size(); i++) {
    if (photons[i]->position.x > maxX) maxX = photons[i]->position.x;
    if (photons[i]->position.x < minX) minX = photons[i]->position.x;

    if (photons[i]->position.y > maxY) maxY = photons[i]->position.y;
    if (photons[i]->position.y < minY) minY = photons[i]->position.y;

    if (photons[i]->position.z > maxZ) maxZ = photons[i]->position.z;
    if (photons[i]->position.z < minZ) minZ = photons[i]->position.z;
  }

  int dimension = -1;

  float rangeX = maxX - minX;
  float rangeY = maxY - minY;
  float rangeZ = maxZ - minZ;

  if (rangeX > rangeY && rangeX > rangeZ) dimension = 0;
  else if (rangeY > rangeX && rangeY > rangeZ) dimension = 1;
  else if (rangeZ > rangeX && rangeZ > rangeY) dimension = 2;
  else dimension = 0;

  return dimension;
}

KDTree* BalanceTree( vector<Photon*> photonPointers ) {

  // Find the dimension where distance is max
  currentMaxDimension = MaxDimension(photonPointers);

  // Sort in the dimension with the maximum distance
  sort(photonPointers.begin(), photonPointers.end(), [](const Photon* lhs, const Photon* rhs) {
    return lhs->position[GetCurrentMaxDimension()] < rhs->position[GetCurrentMaxDimension()];
  });

  KDTree* tree = (KDTree*)malloc(1 * sizeof(KDTree));

  // Find median photon, leftTree and rightTree
  int medianIndex = floor(photonPointers.size() / 2);
  vector<Photon*> leftTree;
  vector<Photon*> rightTree;

  for (int i = 0; i < photonPointers.size(); i++) {
    if (i < medianIndex) leftTree.push_back(photonPointers[i]);
    else if (i > medianIndex) rightTree.push_back(photonPointers[i]);
  }

  // If one more photon left to balance then assign node to it and return
  if (photonPointers.size() == 1) {
    tree->node = photonPointers[0];
    tree->left = endTree;
    tree->right = endTree;

    // tree->left->node = endPhoton;
    // tree->right->node = endPhoton;


    return tree;
  }

  tree->node = photonPointers[medianIndex];

  // Recurse down left and right paths
  if (leftTree.size() != 0) tree->left = BalanceTree(leftTree);
  else tree->node = endPhoton;

  if (rightTree.size() != 0) tree->right = BalanceTree(rightTree);
  else tree->node = endPhoton;

  // Return balance kd tree
  return tree;

}

void Heapify(vector<Photon*>& maxHeap, int heapSize, int i, vec3 position, vec3 normal) {

  int largest = i; // Initialize largest as root node
  int l = 2*i + 1; // left node
  int r = 2*i + 2; // right node

  std::cout << "heapSize: " << heapSize << "or is it? " << maxHeap.size() << '\n';
  std::cout << "before dist" << '\n';

  float distToNode = glm::distance(normal, (position - maxHeap[largest]->position));;
  float distToLeft = glm::distance(normal, (position - maxHeap[l]->position));
  float distToRight = glm::distance(normal, (position - maxHeap[r]->position));

  std::cout << "here" << '\n';

  // If left child is larger than root
  if (l < heapSize && distToLeft > distToNode)
      largest = l;

  // If right child is larger than largest so far
  if (r < heapSize && distToRight > distToNode)
      largest = r;

  // If largest is not root
  if (largest != i)
  {
      std::swap(maxHeap[i], maxHeap[largest]);

      // Recursively heapify the affected sub-tree
      Heapify(maxHeap, heapSize, largest, position, normal);
  }
}

void MaxHeapSort( vector<Photon*>& maxHeap, const Intersection intersection) {
  vec3 position = vec3(intersection.position);
  vec3 normal = vec3(intersection.normal);


  // Build heap (rearrange array)
  for (int i = maxHeap.size() / 2 - 1; i >= 0; i--) {
    std::cout << "heapify build called here with heapSize = " << maxHeap.size() << '\n';
    Heapify(maxHeap, maxHeap.size(), i, position, normal);
  }

  // One by one extract an element from heap
  for (int i = maxHeap.size() - 1; i >= 0; i--)
  {
      std::cout << "heapify extraction called here with heapSize = " << maxHeap.size() << '\n';

      // Move current root to end
      swap(maxHeap[0], maxHeap[i]);

      // call max heapify on the reduced heap
      Heapify(maxHeap, i, 0, position, normal);
  }
}

void LocatePhotons( KDTree* tree, vector<Photon*>& maxHeap, Intersection intersection, float& maxSearchDistance ) {
  float distance;
  vec3 normal = vec3(intersection.normal);
  vec3 position = vec3(intersection.position);
  float squaredDistance = 0.0f;

  // Check signed distance to plane
  if (tree->node != endPhoton) {
    distance = glm::dot(normal, (position - tree->node->position));
    // distance = glm::dot(intersection.normal, (intersection.position - tree->node->position));
  }
  else return;

  // Check if tree has children
  if (tree->left != endTree || tree->right != endTree ) {

    squaredDistance = distance * distance;

    // On left of the plane, search left subtree first
    if (distance < 0) {

      if (tree->left != endTree) {
        LocatePhotons(tree->left, maxHeap, intersection, maxSearchDistance);
      }

      if ((squaredDistance) < (maxSearchDistance)) {
        if (tree->right != endTree) {
          LocatePhotons(tree->right, maxHeap, intersection, maxSearchDistance);
        }
      }
    }
    // On right of the plane, search right subtree first
    else {
      if (tree->right != endTree) {
        LocatePhotons(tree->right, maxHeap, intersection, maxSearchDistance);
      }

      if ((squaredDistance) < (maxSearchDistance)) {

        if (tree->left != endTree) {
          LocatePhotons(tree->left, maxHeap, intersection, maxSearchDistance);
        }
      }
    }
  }

  // Check if distance to this photon is less than photon furthest away
  if ( squaredDistance < maxSearchDistance) {

    // Insert this photon to MaxHeap
    // If heap already has max no of photons, maxHeapify and replace root photon
    if (maxHeap.size() == radianceCount) {
      // Sort MaxHeap to sort photons by distance, root photon being furthest
      MaxHeapSort(maxHeap, intersection);

      maxHeap[0] = tree->node;
    }
    // Else continue adding to the heap
    else {
      maxHeap.push_back(tree->node);

      if (maxHeap.size() > 20) MaxHeapSort(maxHeap, intersection);
    }

    maxSearchDistance = glm::distance(vec3(intersection.position), maxHeap[0]->position);
    maxSearchDistance = maxSearchDistance * maxSearchDistance;
  }
}
