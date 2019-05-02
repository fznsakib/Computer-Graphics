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
#define MAX_RAY_DEPTH 5

// Done
// - Anti aliasing
// - Cramer's rule
// - Load sphere
// - Convert point light to area
// - Reflection

// To do
// - Refraction
// - Build kd tree
// - Photon mapping
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
  vec4 position;
  vec4 direction;
  vec3 power;
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

vector<Triangle> triangles;
vector<Sphere> spheres;

vector<Photon> globalPhotonMap;

vector<Photon*> globalPhotonPointers;

KDTree* globalPhotonTree;

int currentMaxDimension;
int radianceCount;

/* ----------------------------------------------------------------------------*/
/* FUNCTIONS                                                                   */
/* ----------------------------------------------------------------------------*/


bool Update();
void Draw(screen* screen);
bool ClosestIntersection(vec4 start, vec4 dir, Intersection& closestIntersection, int depth );
vec3 DirectLight( const Intersection &i, Light light );
vec3 IndirectLight( const vec3 pixelColour, const vec3 objectColour, const vec3 light );
void PhotonEmission( const int noOfPhotons );
void TracePhoton( Photon& photon );
void GetReflectedDirection(const vec4& incident, const vec4& normal, vec4& reflected);
void GetRefractedDirection(const vec4& incident, const vec4& normal, const float ior, vec4& refracted);
vec4 GenerateRandomDirection();
KDTree* BalanceTree( vector<Photon*> photonPointers );
void LocatePhotons( Photon* rootNode, vector<Photon*> maxHeap, Intersection intersection, float& maxSearchDistance );
int MaxDimension( vector<Photon> photons );
int GetCurrentMaxDimension();
float RandomFloat(float min, float max);
float mix(const float &a, const float &b, const float &mix);
void MaxHeapify( vector<Photon*>& maxHeap, const Intersection intersection);


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

  // KDTree* tree = malloc(1 * sizeof(KDTree));
  KDTree tree;

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
    tree.node = photonPointers[0];
    // return &tree;
  }

  tree.node = photonPointers[medianIndex];

  // Recurse down left and right paths
  if (leftTree.size() != 0) tree.left = BalanceTree(leftTree);
  if (rightTree.size() != 0) tree.right = BalanceTree(rightTree);

  // Return balance kd tree
  return &tree;

}

void MaxHeapify( vector<Photon*>& maxHeap, const Intersection intersection) {
  int swaps = 0;

  do {
    swaps = 0;
    for (int i = floor(maxHeap.size()/2); i >= 0 ; i--) {
      // Find node, leftChild and rightChild's distances to intersection
      int leftIndex, rightIndex;

      if (i == 0) leftIndex = 1;
      else leftIndex = 2 * i;

      if (i == 0) rightIndex = 2;
      else rightIndex = (2 * i) + 1;

      float distToNode  = glm::distance(intersection.position, maxHeap[i]->position);
      float distToLeft  = glm::distance(intersection.position, maxHeap[leftIndex]->position);
      float distToRight = glm::distance(intersection.position, maxHeap[rightIndex]->position);

      Photon* buffer;

      if (distToLeft >= distToRight) {
        if (distToLeft > distToNode) {
          buffer = maxHeap[i];
          maxHeap[i] = maxHeap[leftIndex];
          maxHeap[leftIndex] = buffer;
          swaps += 1;
        }
      }
      else if (distToRight > distToLeft) {
        if (distToRight > distToNode) {
          buffer = maxHeap[i];
          maxHeap[i] = maxHeap[rightIndex];
          maxHeap[rightIndex] = buffer;
          swaps += 1;
        }
      }
    }
  } while(swaps > 0);
}

void LocatePhotons( KDTree* tree, vector<Photon*>& maxHeap, Intersection intersection, float& maxSearchDistance ) {
  // Check if tree has children
  float distance = glm::dot(intersection.normal, (intersection.position - tree->node->position));

  std::cout << "before check" << '\n';
  if (tree->left->node != NULL || tree->right->node != NULL ) {
    std::cout << "after check" << '\n';

    // Check signed distance to plane
    // distance = glm::dot(intersection.normal, (intersection.position - tree->node->position));

    // On left of the plane, search left subtree first
    if (distance < 0 && tree->left->node != NULL) {
      std::cout << "search left subtree (primary)" << '\n';
      LocatePhotons(tree->left, maxHeap, intersection, maxSearchDistance);

      if ((distance * distance) < (maxSearchDistance * maxSearchDistance)) {
        std::cout << "search right subtree (secondary)" << '\n';

        LocatePhotons(tree->right, maxHeap, intersection, maxSearchDistance);
      }
    }
    // On right of the plane, search right subtree first
    else if (distance > 0 && tree->right->node != NULL) {
      std::cout << "search right subtree (primary)" << '\n';

      LocatePhotons(tree->right, maxHeap, intersection, maxSearchDistance);

      if ((distance * distance) < (maxSearchDistance * maxSearchDistance)) {
        std::cout << "search left subtree (secondary)" << '\n';

        LocatePhotons(tree->left, maxHeap, intersection, maxSearchDistance);
      }
    }
  }

  std::cout << "Final node reached" << '\n';

  float squaredDistance = distance * distance;
  // Check if distance to this photon is less than photon furthest away
  if ( squaredDistance < pow(maxSearchDistance, 2)) {

    // If heap already has max no of photons, replace root photon and maxHeapify
    if (maxHeap.size() == radianceCount) {
      maxHeap[0] = tree->node;
    }
    // Else continue adding to the heap
    else {
      maxHeap.push_back(tree->node);
    }

    // Sort MaxHeap to sort photons by distance, root photon being furthest
    MaxHeapify(maxHeap, intersection);

    maxSearchDistance = glm::distance(intersection.position, maxHeap[0]->position);
  }

  std::cout << "Locate complete" << '\n';
}


int main( int argc, char* argv[] )
{

  screen *screen = InitializeSDL( SCREEN_WIDTH, SCREEN_HEIGHT, FULLSCREEN_MODE );

  // Initialise lights
  Light light;
  light.position = vec4( 0, -0.2, -0.7, 1.0 );
  light.colour = vec3( 14.0f * vec3( 1, 1, 1 ));
  lights.push_back(light);

  // Initialise surfaces
  LoadTestModel( triangles, spheres );

  // Initialise photon map before scene rendering
  // int noOfPhotons = 500000;
  int noOfPhotons = 2000;
  PhotonEmission(noOfPhotons);

  // Populate vector with pointers to photons
  for (int i = 0; i < globalPhotonMap.size(); i++) {
    globalPhotonPointers.push_back(&globalPhotonMap[i]);
  }

  // Create and balance photon kd tree
  // globalPhotonTree = (KDTree*)malloc(sizeof(KDTree));
  //
  // globalPhotonTree = BalanceTree(globalPhotonPointers);

  // std::cout << globalPhotonTree << '\n';
  // std::cout << globalPhotonTree->node << '\n';
  // std::cout << globalPhotonTree->node->position.x << '\n';
  // std::cout << globalPhotonTree->left->node->position.x << '\n';

  // std::cout << test->node->position.y << '\n';

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
  // std::cout << "before" << '\n';
  // std::cout << globalPhotonTree->left->node->position.x << '\n';
  // std::cout << (*globalPhotonTree.left->node).direction.x << '\n';
  // std::cout << "after" << '\n';


  /* Clear buffer */
  memset(screen->buffer, 0, screen->height*screen->width*sizeof(uint32_t));

  vec3 black(0.0,0.0,0.0);
  // vec3 pixelColour = vec3( 1.0, 1.0, 1.0 );


  // u and v are coordinates on the 2D screen
  for (int v = 0; v < SCREEN_HEIGHT; v++) {
    for (int u = 0; u < SCREEN_WIDTH; u++) {
      vec4 dir = (vec4((u - SCREEN_WIDTH/2), v - SCREEN_HEIGHT/2, focalLength, 1.0));
      // Calculate the new direction vector after yaw rotation
      dir = R * dir;
      glm::normalize(dir);

      // Apply global photon mapping
      // Intersection intersection;
      // vec3 colour;
      // bool lightCheck = false;
      // vec3 totalRadiance;

      // if (ClosestIntersection(cameraPos, dir, intersection)) {
      //   float distance;
      //   float radius = 0.2f;
      //
      //   // Check if intersecting with light emitter
      //   lightCheck = intersection.position.x >  -0.2f  &&
      //                intersection.position.x <   0.2f  &&
      //                intersection.position.y >= -1.0f  &&
      //                intersection.position.y <  -0.95f &&
      //                intersection.position.z >  -0.8f  &&
      //                intersection.position.z <  -0.4f;


        // Locate photons

        // for (int i = 0; i < globalPhotonMap.size(); i++) {
        //   distance = glm::distance(intersection.position, globalPhotonMap[i].position);
        //
        //   if (distance < radius) colour += globalPhotonMap[i].power;
        // }

        // vector<Photon*> nearestPhotons;
        // float searchRadius = 0.5f;

        // std::cout << "before locating" << '\n';
        // std::cout << globalPhotonTree->node->position.x << '\n';
        // LocatePhotons( globalPhotonTree, nearestPhotons, intersection, searchRadius );

        // Calculate radiance estimate
      //   for (int i = 0; i < nearestPhotons.size(); i++ ) {
      //     totalRadiance += nearestPhotons[i]->power;
      //   }
      //
      //   totalRadiance = totalRadiance / float((M_PI * (searchRadius * searchRadius)));
      // }
      // if (lightCheck) {
      //   PutPixelSDL(screen, u, v, vec3(1.0f, 1.0f, 1.0f));
      // }
      // else PutPixelSDL(screen, u, v, totalRadiance);


      // OLD RENDERING METHOD
      vec3 pixelColour = vec3(0.0f, 0.0f, 0.0f);
      vec3 indirectLight = 0.5f * vec3(1, 1, 1);
      bool validRay = false;

      // Trace multiple rays around each ray and average color for anti aliasing
      // Change to 4 x 4 later
      for (int i = -2; i <= 2; ++i) {
        for (int j = -2; j <= 2; ++j ) {
          // Multiply i and j by distance to surrounding rays being used
          float multiplier = 0.25;
          vec4 newDir = vec4(dir.x + (multiplier * i), dir.y + (multiplier * j), focalLength, 1.0);

          Intersection intersection;
          int depth = 0;
          bool closestIntersection = ClosestIntersection(cameraPos, newDir, intersection, depth);

          // If light intersects with object then draw it
          if (closestIntersection == true) {
            validRay = true;
            vec3 objectColour = intersection.colour;

            // Direct illumination
            for (int i = 0; i < lights.size(); i++) {
              if (DirectLight( intersection, lights[i] ) != vec3(0.0f, 0.0f, 0.0f)) {
                pixelColour += DirectLight( intersection, lights[i] );
              }
            }

            // Indirect illumination for diffuse objects
            if (intersection.material[1] == 0.0f) {
              pixelColour = IndirectLight(pixelColour, objectColour, indirectLight);
            }
          }
        }
      }
      if (validRay == true) {
        // Average light for the multiple rays for anti aliasing
        vec3 averageLight = pixelColour/16.0f;
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


bool ClosestIntersection( vec4 start, vec4 dir,
                          Intersection& closestIntersection, int depth) {

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

    // Reflection
    // if (closestIntersection.material[1] > 0.0f && closestIntersection.material[2] == 0.0f) {
    //   Intersection reflection;
    //
    //   // Calculate new reflected direction of ray
    //   vec4 reflectionDir = glm::reflect(dir, closestIntersection.normal);
    //   vec3 indirectLight = 0.5f * vec3(1, 1, 1);
    //
    //   // Trace reflected ray
    //   if (ClosestIntersection(closestIntersection.position + (closestIntersection.normal * 0.00001f), reflectionDir, reflection, 0)) {
    //     vec3 colour = IndirectLight(closestIntersection.colour, reflection.colour, indirectLight);
    //     reflection.colour = colour;
    //   }
    //   else {
    //     return false;
    //   }
    //
    //   closestIntersection.colour = reflection.colour;
    //
    //   return true;
    // }


    // REFLECTION AND REFRACTION
    // If the normal and the view direction are not opposite to each other
    // reverse the normal direction. That also means we are inside the sphere so set
    // the inside bool to true. Finally reverse the sign of IdotN which we want
    // positive.
    bool inside = false;
    if (glm::dot(dir, closestIntersection.normal) > 0) {
      closestIntersection.normal = -closestIntersection.normal;
      inside = true;
    }

    vec3 indirectLight = 0.5f * vec3(1, 1, 1);
    float specular = closestIntersection.material[1];
    float transmission = closestIntersection.material[2];
    Intersection reflection;
    Intersection refraction;

    // Check if reflective or refractive
    if ((specular > 0.0f || transmission > 0.0f) && depth < MAX_RAY_DEPTH) {
      float facingratio = -glm::dot(dir, closestIntersection.normal);

      // Reflection/refraction mix value
      float fresnelEffect = mix(pow(1 - facingratio, 3), 1, 0.1);

      // Calculate new reflected direction of ray
      vec4 reflectionDir;
      GetReflectedDirection(dir, closestIntersection.normal, reflectionDir);

      // Trace reflected ray
      depth += 1;
      if (ClosestIntersection(closestIntersection.position + (closestIntersection.normal * 0.00001f), reflectionDir, reflection, depth)) {
        vec3 colour = IndirectLight(closestIntersection.colour, reflection.colour, indirectLight);
        reflection.colour = colour;
      }

      // Check if refractive
      if (transmission > 0.0f) {
        // Compute index of refraction depending on whether light is inside or outside object
        float ior = 1.5;
        float eta = (inside) ? ior : 1 / ior;

        // Calculate new refracted direction of ray
        vec4 refractionDir;
        GetRefractedDirection(dir, closestIntersection.normal, eta, refractionDir);
        glm::normalize(refractionDir);

        // Trace refracted ray
        depth += 1;
        ClosestIntersection(closestIntersection.position - (closestIntersection.normal * 0.00001f), refractionDir, refraction, depth);
      }

      // std::cout << "reflection colour: " << reflection.colour.x << reflection.colour.y << reflection.colour.z << '\n';
      // std::cout << "refraction colour: " << refraction.colour.x << refraction.colour.y << refraction.colour.z << '\n';

      // std::cout << fresnelEffect << '\n';

      // Now compute the mix of reflection and refraction colours
      // vec3 computedColour = ((reflection.colour * fresnelEffect) +
      //                       (refraction.colour * (1 - fresnelEffect) * transmission)
      //                       * closestIntersection.colour);

      // std::cout << refraction.colour.x << '\n';
      vec3 computedColour = reflection.colour + refraction.colour;
      // vec3 computedColour = refraction.colour;

      // std::cout << computedColour.x << computedColour.y << computedColour.z << '\n';

      closestIntersection.colour += computedColour;
    }

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

  // Check if distance less than to light. Emit ray from a small distance off the surface
  // so that it doesn't intersect with itself
  if (ClosestIntersection(i.position + (normal * 0.00001f), direction, intersection, 0)) {
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

    // Rejection sampling to find random photon direction
    do {
      x = RandomFloat(-1.0, 1.0);
      y = RandomFloat(-1.0, 1.0);
      z = RandomFloat(-1.0, 1.0);
    } while ((x*x) + (y*y) + (z*z) > 1);

    // Change to square light later
    photon.position = vec4(RandomFloat(-0.2f, 0.2f), -0.95f, RandomFloat(-0.8f, -0.4f), 1.0f);
    float u = RandomFloat(0.0f, 1.0f);
    float v = 2 * M_PI * RandomFloat(0.0f, 1.0f);
    photon.direction = vec4((cos(v) * sqrt(u)), sqrt(1- u), (sin(v) * sqrt(u)), 1.0f);

    // photon.position = vec4(lights[0].position.x, lights[0].position.y, lights[0].position.z, 1.0f);
    // photon.direction = vec4(x, y, z, 1.0);
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

    if (!ClosestIntersection(photon.position, photon.direction, intersection, 0)) {
      break;
    }

    // Use Russian Roulette on material properties to decide to photon action
    float random = RandomFloat(0.0f, 1.0f);

    // Add photon for each intersection
    // DIFFUSE REFLECTION
    if (random < intersection.material[0]) {
      // Reflect photon with new position in new direction
      photon.position = intersection.position;
      photon.power = photon.power * intersection.colour;

      globalPhotonMap.push_back(photon);

      // Reflect in random direction
      vec4 randomDirection = GenerateRandomDirection();
      Intersection intersectionTest;

      while (!ClosestIntersection(photon.position, randomDirection, intersectionTest, 0)) {
        randomDirection = GenerateRandomDirection();
      }

      photon.direction = randomDirection;
    }
    // SPECULAR REFLECTION
    else if (random < (intersection.material[0] + intersection.material[1])) {
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


float mix(const float &a, const float &b, const float &mix) {
    return b * mix + a * (1 - mix);
}
