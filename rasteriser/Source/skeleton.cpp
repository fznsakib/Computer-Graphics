#include <iostream>
#include <algorithm>
#include <glm/glm.hpp>
#include <SDL.h>
#include "SDLauxiliary.h"
#include "TestModelH.h"
#include <stdint.h>

using namespace std;
using glm::vec3;
using glm::mat3;
using glm::vec4;
using glm::mat4;

SDL_Event event;

#define SCREEN_WIDTH 640
#define SCREEN_HEIGHT 512
#define FULLSCREEN_MODE false


/* ----------------------------------------------------------------------------*/
/* VARIABLES                                                                   */
/* ----------------------------------------------------------------------------*/

float focalLength = 512;
vec4 cameraPos( 0, 0, -3.001, 1 );
// float yaw = 0.0;
// mat4 R(1.0f);
vec3 currentColor;
// Stores depth information of screen pixels. This is stored as 1/z.
float depthBuffer[SCREEN_HEIGHT][SCREEN_WIDTH];
vec3 screenBuffer[SCREEN_HEIGHT][SCREEN_WIDTH];
int shadowBuffer[SCREEN_HEIGHT][SCREEN_WIDTH];


// Initialise light variables
vec4 lightPos(0, -0.5, 0, 1);
vec3 lightPower = 11.0f*vec3( 1, 1, 1 );
vec3 indirectLightPowerPerArea = 0.08f*vec3( 1, 1, 1 );

// Normal and reflectance values to be passed to PixelShader
vec4 currentNormal;
vec3 currentReflectance;

// Choose random colour or correct colour
int randColourSelect = 0;

// Store frames per second to print.
int fps;

struct Pixel
{
  int x;
  int y;
  float zinv;
  //vec3 illumination;
  vec4 pos3d;
};

struct Vertex
{
  vec4 position;
  //vec4 normal;
  //glm::vec3 reflectance;
};

/* ----------------------------------------------------------------------------*/
/* FUNCTIONS                                                                   */
/* ----------------------------------------------------------------------------*/

bool Update();
void Draw(screen* screen);
void VertexShader( const Vertex& v, Pixel& p );
void Interpolate( Pixel a, Pixel b, vector<Pixel>& result );
void DrawLineSDL( screen* screen, vector<glm::ivec2> line, vec3 color );
void DrawPolygon( screen* screen, const vector<Vertex>& vertices);
void ComputePolygonRows( const vector<Pixel>& vertexPixels, vector<Pixel>& leftPixels, vector<Pixel>& rightPixels);
void DrawPolygonRows( screen* screen, const vector<Pixel>& leftPixels, const vector<Pixel>& rightPixels );
void PixelShader( const Pixel& p, screen* screen );
vec3 calculateIllumination( const Pixel& p, vec4 currentNormal );
void toClipSpace( vector<Triangle>& v );
vector<Triangle> clipZ( vector<Triangle> triangles );
vector<Triangle> clip( vector<Triangle>& triangles, int plane );
vector<Triangle> createShadowVolume( vector<Triangle>& triangles );
void toCameraSpace( vector<Triangle>& triangles );
float surroundingShadowSum( int x, int y );

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

  // /* TO TEST INTERPOLATION VALUES */
  // vector<Pixel> vertexPixels(3);
  // vertexPixels[0] = (Pixel) {.x = 10, .y = 5};
  // vertexPixels[1] = (Pixel) {.x = 5, .y = 10};
  // vertexPixels[2] = (Pixel) {.x = 15, .y = 15};
  // vector<Pixel> leftPixels;
  // vector<Pixel> rightPixels;
  // ComputePolygonRows( vertexPixels, leftPixels, rightPixels );
  // for( int row=0; row<leftPixels.size(); ++row )
  // {
  //   cout << "Start: ("
  //   << leftPixels[row].x << ","
  //   << leftPixels[row].y << "). "
  //   << "End: ("
  //   << rightPixels[row].x << ","
  //   << rightPixels[row].y << "). " << endl;
  // }
}

/*Place your drawing here*/
void Draw(screen* screen)
{
  vector<Triangle> triangles;
  vector<Triangle> boxes;
  LoadTestModel( triangles, boxes );
  //std::vector<glm::ivec2> points;

  // Transform triangles to camera space.
  toCameraSpace(triangles);
  toCameraSpace(boxes);

  // Take all subjects in the room and add shadows.
   boxes = createShadowVolume(boxes);

   // Add the shadowed boxes triangles to the main triangles vector.
   for (int i = 0; i < boxes.size(); ++i) {
     triangles.push_back(boxes[i]);
   }

  // Transform triangles to clip space
  toClipSpace(triangles);

  // Clip the list of triangles on all planes of the frustrum
  // 1: Left of camera  2: Right of camera  3: Bottom of camera  4: Top of camera
  // 5: Behind camera
  vector<Triangle> clippedTriangles = clip(triangles, 1);
  clippedTriangles = clip(clippedTriangles, 2);
  clippedTriangles = clip(clippedTriangles, 3);
  clippedTriangles = clip(clippedTriangles, 4);
  clippedTriangles = clip(clippedTriangles, 5);
  clippedTriangles = clip(clippedTriangles, 6);
  printf("%lu\n",clippedTriangles.size() );

  /* Clear buffer */
  memset(screen->buffer, 0, screen->height*screen->width*sizeof(uint32_t));

  // Clear depth buffer.
  memset(depthBuffer, 0, SCREEN_WIDTH * SCREEN_HEIGHT * sizeof(float));

  // Clear screen buffer.
  memset(screenBuffer, 0, SCREEN_WIDTH * SCREEN_HEIGHT * sizeof(vec3));

  // Clear shadow buffer.
  memset(shadowBuffer, 0, SCREEN_WIDTH * SCREEN_HEIGHT * sizeof(int));

  // Go through all clipped triangles
  for( uint32_t i=0; i<clippedTriangles.size(); ++i ) {

     // Load triangle information into vertices
     vector<Vertex> vertices(3);

     vertices[0].position = clippedTriangles[i].v0;
     vertices[1].position = clippedTriangles[i].v1;
     vertices[2].position = clippedTriangles[i].v2;

     // Set current colour, reflectance and normal to the colour of the triangle.
     currentColor = clippedTriangles[i].color;
     currentNormal = clippedTriangles[i].normal;
     //currentReflectance = ??;

     // Draw and colour the triangles.
     DrawPolygon( screen, vertices );
  }
  // Go through the screen buffer and shadow buffer, and decide what to print as pixels/shadows.
  for (int y = 1; y < SCREEN_HEIGHT-1; ++y) {
    for (int x = 1; x < SCREEN_WIDTH-1; ++x) {
      // Go through shadow buffer first and reduce illumination at shadow areas.
      if (shadowBuffer[y][x] == 1) {
        // Soft shadows first.
        if (surroundingShadowSum(y,x) < 0.6) {
          screenBuffer[y][x] -= vec3(0.05f, 0.05f, 0.05f);
        }
        else if (surroundingShadowSum(y,x) < 0.7) {
          screenBuffer[y][x] -= vec3(0.08f, 0.08f, 0.08f);
        }
         else if (surroundingShadowSum(y,x) < 0.8) {
          screenBuffer[y][x] -= vec3(0.1f, 0.1f, 0.1f);
        }
        else if (surroundingShadowSum(y,x) < 0.9) {
         screenBuffer[y][x] -= vec3(0.12f, 0.12f, 0.12f);
       }
       else {
         screenBuffer[y][x] -= vec3(0.5f, 0.5f, 0.5f);
       }
      }
      PutPixelSDL(screen, x, y, (screenBuffer[y][x] +
        screenBuffer[y-1][x] + screenBuffer[y+1][x] + screenBuffer[y][x-1] +
        screenBuffer[y][x+1] + screenBuffer[y-1][x-1] + screenBuffer[y+1][x+1] +
        screenBuffer[y+1][x-1] + screenBuffer[y-1][x+1])/9.0f);
    }
  }
}

/*Place updates of parameters here*/
bool Update()
{
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
        // LIGHT MOVEMENT
        case SDLK_w:
          lightPos += vec4(0, 0, 0.1, 0);
        break;
        case SDLK_s:
          lightPos += vec4(0, 0, -0.1, 0);
        break;
        case SDLK_a:
          lightPos += vec4(-0.1, 0, 0, 0);
        break;
        case SDLK_d:
          lightPos += vec4(0.1, 0, 0, 0);
        break;
        case SDLK_q:
          lightPos += vec4(0, -0.1, 0, 0);
        break;
        case SDLK_e:
          lightPos += vec4(0, 0.1, 0, 0);
        break;
        /* Reduce surrounding brightness */
        case SDLK_1:
        printf("YES");
          indirectLightPowerPerArea -= vec3(0.005f, 0.005f, 0.005f);
        break;
        /* Increase surrounding brightness */
        case SDLK_2:
          indirectLightPowerPerArea += vec3(0.005f, 0.005f, 0.005f);
        break;
        // CAMERA MOVEMENT
        case SDLK_UP:
        /* Move camera forward */
        cameraPos += vec4(0, 0, 0.1, 0);
        break;
        case SDLK_DOWN:
        /* Move camera backwards */
        cameraPos += vec4(0, 0, -0.1, 0);
        break;
        case SDLK_LEFT:
        /* Move camera left */
        cameraPos += vec4(-0.1, 0, 0, 0);
        break;
        case SDLK_RIGHT:
        /* Move camera right */
        cameraPos += vec4(0.1, 0, 0, 0);
        break;
        case SDLK_z:
        /* Move camera down */
        cameraPos += vec4(0, -0.1, 0, 0);
        break;
        case SDLK_x:
        /* Move camera up */
        cameraPos += vec4(0, 0.1, 0, 0);
        break;
        case SDLK_f:
        /* Increase focal length */
        focalLength +=5;
        break;
        case SDLK_g:
        /* Decrease focal length */
        focalLength -=5;
        break;
        case SDLK_SPACE:
        /* Toggle colour system */
        randColourSelect = (randColourSelect+1)%3;
        break;
        case SDLK_ESCAPE:
        /* Move camera quit */
        return false;
      }
    }
  }
  return true;
}

/*----------------- ADDED FUNCTION TO SHADE IN POLYGONS --------------------- */
void DrawPolygon(screen* screen, const vector<Vertex>& vertices) {
  int V = vertices.size();

  vector<Pixel> vertexPixels(V);
  for (int i=0; i<V; ++i) {
    VertexShader( vertices[i], vertexPixels[i] );
  }
  vector<Pixel> leftPixels;
  vector<Pixel> rightPixels;
  ComputePolygonRows( vertexPixels, leftPixels, rightPixels );
  DrawPolygonRows( screen, leftPixels, rightPixels );
}

void ComputePolygonRows( const vector<Pixel>& vertexPixels, vector<Pixel>& leftPixels, vector<Pixel>& rightPixels) {
  int max = -numeric_limits<int>::max();
  int min = +numeric_limits<int>::max();
  int rows = 0;

  // Find min and max y-value of polygon and compute number of rows it occupies.
  for (int i = 0; i < vertexPixels.size(); ++i) {
    Pixel currentVertex = vertexPixels[i];
    if (currentVertex.y > max) {
      max = currentVertex.y;
    }
    if (currentVertex.y < min) {
      min = currentVertex.y;
    }
  }
  rows = (max - min) + 1;

  // Resize leftPixels and rightPixels to fit the correct amount of rows
  leftPixels.resize(rows);
  rightPixels.resize(rows);

  // Initialize x coordinates in leftPixels to large vals and rightPixels to small vals.
  for (int j = 0; j < rows; ++j) {
    leftPixels[j].x = +numeric_limits<int>::max();
    leftPixels[j].y = min + j;
    rightPixels[j].x = -numeric_limits<int>::max();
    rightPixels[j].y = min + j;
  }

  //Loop through all edges, use linear interpolation to find x-coordinate for
  //each row it occupies. Update the corresponding vals in rightPixels and leftPixels
  Pixel a;
  Pixel b;
  for(int i = 0; i < vertexPixels.size(); i++) {
    a = vertexPixels[i];
    if (i == vertexPixels.size() - 1) b = vertexPixels[0];
    else b = vertexPixels[i + 1];

    // Calculate number of pixels to draw in interpolation
    Pixel delta;
    delta.x = glm::abs( a.x - b.x );
    delta.y = glm::abs( a.y - b.y );
    int pixels = glm::max( delta.x, delta.y ) + 1;

    // Populate vector with pixels to be drawn
    vector<Pixel> line( pixels );
    Interpolate( a, b, line );

    /* interpolation changed to floor() instead of round() to achieve correct vals. */
    /*----------------------------------------------*/
    for (int j = 0; j < line.size(); ++j) {
      int inTriangleHeightIndex = line[j].y - min;
      if (line[j].x <= leftPixels[ inTriangleHeightIndex ].x && line[j].y - min >= 0) { // SEG FAULT FIXED HERE
        leftPixels[ inTriangleHeightIndex ].x = line[j].x;
        leftPixels[ inTriangleHeightIndex ].zinv = line[j].zinv;
        leftPixels[ inTriangleHeightIndex ].pos3d = line[j].pos3d;
      }
      if (line[j].x >= rightPixels[ inTriangleHeightIndex ].x && line[j].y - min >= 0) { // SEG FAULT FIXED HERE
        rightPixels[ inTriangleHeightIndex ].x = line[j].x;
        rightPixels[ inTriangleHeightIndex ].zinv = line[j].zinv;
        rightPixels[ inTriangleHeightIndex ].pos3d = line[j].pos3d;
      }
    }
    /*----------------------------------------------*/
  }
}

void DrawPolygonRows( screen* screen, const vector<Pixel>& leftPixels, const vector<Pixel>& rightPixels ) {
  for (int y = 0; y < leftPixels.size(); ++y) {
    vector<Pixel> currentLine(rightPixels[y].x - leftPixels[y].x + 1);
    Interpolate(leftPixels[y], rightPixels[y], currentLine);
    for (int x = 0; x < (rightPixels[y].x - leftPixels[y].x); ++x) {
      PixelShader(currentLine[x], screen );
    }
  }
}

void VertexShader( const Vertex& v, Pixel& p ) {

  float x = (focalLength * (v.position.x / v.position.z)) + (SCREEN_WIDTH/2);
  float y = (focalLength * (v.position.y / v.position.z)) + (SCREEN_HEIGHT/2);

  int x2D = static_cast<int>(x);
  int y2D = static_cast<int>(y);

  p.x = x2D;
  p.y = y2D;
  p.zinv = 1/v.position.z;

  // Compute illumination (vertex illumination)
  // vec4 r = lightPos - v.position;
  // vec3 vec3r = vec3(r);
  // float r_magnitude = pow(r[0],2) + pow(r[1],2) + pow(r[2],2);
  //
  // vec3 normal = vec3(v.normal);
  //
  // float vecProduct = glm::dot(vec3r, normal);
  // vec3 D = (lightPower * glm::max(vecProduct, 0.0f)) / (float)(4.0f * M_PI * r_magnitude);
  //
  // vec3 totalIllumination = D + indirectLightPowerPerArea;
  //
  // p.illumination = totalIllumination;

  // Store 3D position to pixel for illumination (per pixel illumination)
  p.pos3d = v.position;
}

void Interpolate( Pixel a, Pixel b, vector<Pixel>& result ) {
  // Perspective correction
  a.pos3d.x = a.pos3d.x * a.zinv;
  a.pos3d.y = a.pos3d.y * a.zinv;

  b.pos3d.x = b.pos3d.x * b.zinv;
  b.pos3d.y = b.pos3d.y * b.zinv;
  int N = result.size();

  float step_x = (float)(b.x-a.x) / float(max(N-1,1));
  float step_y = (float)(b.y-a.y) / float(max(N-1,1));
  float step_z = (b.zinv-a.zinv) / float(max(N-1,1));

  float stepPos3dX = (b.pos3d.x - a.pos3d.x) / float(max(N-1,1));
  float stepPos3dY = (b.pos3d.y - a.pos3d.y) / float(max(N-1,1));
  float stepPos3dZ = (b.pos3d.z - a.pos3d.z) / float(max(N-1,1));

  for( int i=0; i<N; ++i ) {
     result[i].x = floor(a.x + (step_x * i));
     result[i].y = floor(a.y + (step_y * i));
     result[i].zinv = (a.zinv + (step_z * i));

     // Divide by z-inverse for perspective correction
     result[i].pos3d.z = 1/result[i].zinv;
     result[i].pos3d.x = (a.pos3d.x + (stepPos3dX * i)) / result[i].zinv;
     result[i].pos3d.y = (a.pos3d.y + (stepPos3dY * i)) / result[i].zinv;
     result[i].pos3d.w = 1.0f;
  }
}

void DrawLineSDL( screen* screen, vector<glm::ivec2> line, vec3 color ) {
    for(int i = 0; i < line.size(); i++) {
      PutPixelSDL(screen, line[i].x, line[i].y, color);
    }
}

void PixelShader( const Pixel& p, screen* screen ) {
  int x = p.x;
  int y = p.y;

  // Random colour generator
  float LO = 0.2f;
  float HI = 0.5f;
  float r0 = LO + static_cast <float> (rand()) / (static_cast <float> (RAND_MAX/HI-LO));
  float r1 = LO + static_cast <float> (rand()) / (static_cast <float> (RAND_MAX/HI-LO));
  float r2 = LO + static_cast <float> (rand()) / (static_cast <float> (RAND_MAX/HI-LO));
  vec3 randColour = vec3(r0, r1, r2);
  vec3 nightVision = vec3(r0-0.2f, 1.0f, r2-0.2f);

  // Draw on screen if it is closer to the camera, and choose colour.
  if (x >= 0 && x < SCREEN_WIDTH && y >= 0 && y < SCREEN_HEIGHT) {
    if (p.zinv >= depthBuffer[y][x] && currentColor.x >= 0) {
      depthBuffer[y][x] = p.zinv;
      switch (randColourSelect) {
        case 0:
        // Choose normal colour
        // PutPixelSDL( screen, x, y, currentColor * calculateIllumination(p, currentNormal) );
        screenBuffer[y][x] = currentColor * calculateIllumination(p, currentNormal);
        break;
        case 1:
        // Choose random colour
        // PutPixelSDL( screen, x, y, randColour * calculateIllumination(p, currentNormal) );
        screenBuffer[y][x] = randColour * calculateIllumination(p, currentNormal);
        break;
        case 2:
        // Choose nightVision
        // PutPixelSDL( screen, x, y, nightVision * calculateIllumination(p, currentNormal) );
        screenBuffer[y][x] = nightVision * calculateIllumination(p, currentNormal);
        break;
      }
    }
    else if (p.zinv > depthBuffer[y][x] && currentColor.x < 0) {
      shadowBuffer[y][x] = 1;
    }
  }
}

vec3 calculateIllumination(const Pixel& p, vec4 currentNormal) {
  vec4 r = (lightPos - cameraPos) - p.pos3d;
  vec3 vec3r = vec3(r);
  float r_magnitude = pow(r[0],2) + pow(r[1],2) + pow(r[2],2);

  vec3 normal = vec3(currentNormal);

  float vecProduct = glm::dot(vec3r, normal);
  vec3 D = (lightPower * glm::max(vecProduct, 0.0f)) / (float)(4.0f * M_PI * r_magnitude);

  vec3 totalIllumination = D + indirectLightPowerPerArea;

  return totalIllumination;

}

/* THIS FUNCTION CONVERTS ALL TRIANGLE INFORMATION FROM 3D TO CLIP SPACE */
void toClipSpace(vector<Triangle>& triangles) {
  for (int i = 0; i<triangles.size(); ++i) {

    // W = Z/f
    triangles[i].v0.w = triangles[i].v0.z/focalLength;
    triangles[i].v1.w = triangles[i].v1.z/focalLength;
    triangles[i].v2.w = triangles[i].v2.z/focalLength;
  }
}

void toCameraSpace(vector<Triangle>& triangles) {
  for (int i = 0; i<triangles.size(); ++i) {
    // Changing from world space to camera space.
    triangles[i].v0 = triangles[i].v0 - cameraPos;
    triangles[i].v0.w = 1.0f;
    triangles[i].v1 = triangles[i].v1 - cameraPos;
    triangles[i].v1.w = 1.0f;
    triangles[i].v2 = triangles[i].v2 - cameraPos;
    triangles[i].v2.w = 1.0f;
  }
}

/* THIS FUNCTION CLIPS THE VIEWING FRUSTRUM. EACH CASE CORRESPONDS TO A DIFFERENT
   PLANE */
vector<Triangle> clip(vector<Triangle>& triangles, int plane) {

  // Vector to store triangles kept after clipping process.
  vector<Triangle> clippedTriangles;

  // Array to store the w-dimension corrected limits for each vertex.
  float dot[3];

  switch (plane) {
    // Clip at left side x-plane
    case 1:
    for (int i = 0; i < triangles.size(); ++i) {
      dot[0] = triangles[i].v0.w * -SCREEN_WIDTH/2; // the line of x = -SCREEN_WIDTH/2 in the clip space
      dot[1] = triangles[i].v1.w * -SCREEN_WIDTH/2;
      dot[2] = triangles[i].v2.w * -SCREEN_WIDTH/2;

      // CASE ALL VERTICES ARE IN THE PLANE
      if (triangles[i].v0.x > dot[0] && triangles[i].v1.x > dot[1] && triangles[i].v2.x > dot[2]) {
        clippedTriangles.push_back(triangles[i]);
      }

      // CASE V0 IS IN
      else if (triangles[i].v0.x > dot[0] && triangles[i].v1.x <= dot[1] && triangles[i].v2.x <= dot[2]) {
        // Get all x and w values.
        float x0 = triangles[i].v0.x;
        float x1 = triangles[i].v1.x;
        float x2 = triangles[i].v2.x;

        float w0 = triangles[i].v0.w;
        float w1 = triangles[i].v1.w;
        float w2 = triangles[i].v2.w;

        // Calculate the line intersection with the plane to change vertex positions within the plane.
        float t_01 = (x0 + SCREEN_WIDTH/2*w0)/(-SCREEN_WIDTH/2*w1 + SCREEN_WIDTH/2*w0 - x1 + x0);
        float t_02 = (x0 + SCREEN_WIDTH/2*w0)/(-SCREEN_WIDTH/2*w2 + SCREEN_WIDTH/2*w0 - x2 + x0);

        // Set new vertex positions.
        triangles[i].v1 = triangles[i].v0 + t_01*(triangles[i].v1-triangles[i].v0);
        triangles[i].v2 = triangles[i].v0 + t_02*(triangles[i].v2-triangles[i].v0);

        // Add new triangle to vector.
        clippedTriangles.push_back(triangles[i]);
      }

      // CASE V1 IS IN
      else if (triangles[i].v0.x <= dot[0] && triangles[i].v1.x > dot[1] && triangles[i].v2.x <= dot[2]) {
        // Get all x and w values.
        float x0 = triangles[i].v0.x;
        float x1 = triangles[i].v1.x;
        float x2 = triangles[i].v2.x;

        float w0 = triangles[i].v0.w;
        float w1 = triangles[i].v1.w;
        float w2 = triangles[i].v2.w;

        // Calculate the line intersection with the plane to change vertex positions within the plane.
        float t_10 = (x1 + SCREEN_WIDTH/2*w1)/(-SCREEN_WIDTH/2*w0 + SCREEN_WIDTH/2*w1 - x0 + x1);
        float t_12 = (x1 + SCREEN_WIDTH/2*w1)/(-SCREEN_WIDTH/2*w2 + SCREEN_WIDTH/2*w1 - x2 + x1);

        // Set new vertex positions.
        triangles[i].v0 = triangles[i].v1 + t_10*(triangles[i].v0-triangles[i].v1);
        triangles[i].v2 = triangles[i].v1 + t_12*(triangles[i].v2-triangles[i].v1);

        // Add new triangle to vector.
        clippedTriangles.push_back(triangles[i]);
      }

      // CASE V2 IS IN
      else if (triangles[i].v0.x <= dot[0] && triangles[i].v1.x <= dot[1] && triangles[i].v2.x > dot[2]) {
        // Get all x and w values.
        float x0 = triangles[i].v0.x;
        float x1 = triangles[i].v1.x;
        float x2 = triangles[i].v2.x;

        float w0 = triangles[i].v0.w;
        float w1 = triangles[i].v1.w;
        float w2 = triangles[i].v2.w;

        // Calculate the line intersection with the plane to change vertex positions within the plane.
        float t_21 = (x2 + SCREEN_WIDTH/2*w2)/(-SCREEN_WIDTH/2*w1 + SCREEN_WIDTH/2*w2 - x1 + x2);
        float t_20 = (x2 + SCREEN_WIDTH/2*w2)/(-SCREEN_WIDTH/2*w0 + SCREEN_WIDTH/2*w2 - x0 + x2);

        // Set new vertex positions.
        triangles[i].v1 = triangles[i].v2 + t_21*(triangles[i].v1-triangles[i].v2);
        triangles[i].v0 = triangles[i].v2 + t_20*(triangles[i].v0-triangles[i].v2);

        // Add new triangle to vector.
        clippedTriangles.push_back(triangles[i]);
      }

              /* WHEN 2 VERTICES ARE INSIDE, SPLIT INTO 2 TRIANGLES */

      // CASE V0 AND V1 ARE IN
      else if (triangles[i].v0.x > dot[0] && triangles[i].v1.x > dot[1] && triangles[i].v2.x <= dot[2]) {
        // Get all x and w values.
        float x0 = triangles[i].v0.x;
        float x1 = triangles[i].v1.x;
        float x2 = triangles[i].v2.x;

        float w0 = triangles[i].v0.w;
        float w1 = triangles[i].v1.w;
        float w2 = triangles[i].v2.w;

        // Calculate the line intersection with the plane to change vertex positions within the plane.
        float t_12 = (x1 + SCREEN_WIDTH/2*w1)/(-SCREEN_WIDTH/2*w2 + SCREEN_WIDTH/2*w1 - x2 + x1);
        float t_02 = (x0 + SCREEN_WIDTH/2*w0)/(-SCREEN_WIDTH/2*w2 + SCREEN_WIDTH/2*w0 - x2 + x0);

        // Store new calculated vertices here.
        vec4 newPoint_12;
        vec4 newPoint_02;
        newPoint_12 = triangles[i].v1 + t_12*(triangles[i].v2-triangles[i].v1);
        newPoint_02 = triangles[i].v0 + t_02*(triangles[i].v2-triangles[i].v0);

        // Set outside vertex to one of the new vertices.
        triangles[i].v2 = newPoint_02;

        // Create a new triangle with the new vertices and one of the inside vertex.
        // Set the normals to be the same.
        Triangle extraTriangle(newPoint_02, newPoint_12, triangles[i].v1, triangles[i].color);
        extraTriangle.normal = triangles[i].normal;

        clippedTriangles.push_back(triangles[i]);
        clippedTriangles.push_back(extraTriangle);
      }

      // CASE V0 AND V2 ARE IN
      else if (triangles[i].v0.x > dot[0] && triangles[i].v1.x <= dot[1] && triangles[i].v2.x > dot[2]) {
        // Get all x and w values.
        float x0 = triangles[i].v0.x;
        float x1 = triangles[i].v1.x;
        float x2 = triangles[i].v2.x;

        float w0 = triangles[i].v0.w;
        float w1 = triangles[i].v1.w;
        float w2 = triangles[i].v2.w;

        // Calculate the line intersection with the plane to change vertex positions within the plane.
        float t_01 = (x0 + SCREEN_WIDTH/2*w0)/(-SCREEN_WIDTH/2*w1 + SCREEN_WIDTH/2*w0 - x1 + x0);
        float t_21 = (x2 + SCREEN_WIDTH/2*w2)/(-SCREEN_WIDTH/2*w1 + SCREEN_WIDTH/2*w2 - x1 + x2);

        // Store new calculated vertices here.
        vec4 newPoint_01;
        vec4 newPoint_21;
        newPoint_01 = triangles[i].v0 + t_01*(triangles[i].v1-triangles[i].v0);
        newPoint_21 = triangles[i].v2 + t_21*(triangles[i].v1-triangles[i].v2);

        // Set outside vertex to one of the new vertices.
        triangles[i].v1 = newPoint_01;

        // Create a new triangle with the new vertices and one of the inside vertex.
        // Set the normals to be the same.
        Triangle extraTriangle(newPoint_01, newPoint_21, triangles[i].v2, triangles[i].color);
        extraTriangle.normal = triangles[i].normal;

        clippedTriangles.push_back(triangles[i]);
        clippedTriangles.push_back(extraTriangle);
      }

      // CASE V1 AND V2 ARE IN
      else if (triangles[i].v0.x <= dot[0] && triangles[i].v1.x > dot[1] && triangles[i].v2.x > dot[2]) {
        // Get all x and w values.
        float x0 = triangles[i].v0.x;
        float x1 = triangles[i].v1.x;
        float x2 = triangles[i].v2.x;

        float w0 = triangles[i].v0.w;
        float w1 = triangles[i].v1.w;
        float w2 = triangles[i].v2.w;

        // Calculate the line intersection with the plane to change vertex positions within the plane.
        float t_10 = (x1 + SCREEN_WIDTH/2*w1)/(-SCREEN_WIDTH/2*w0 + SCREEN_WIDTH/2*w1 - x0 + x1);
        float t_20 = (x2 + SCREEN_WIDTH/2*w2)/(-SCREEN_WIDTH/2*w0 + SCREEN_WIDTH/2*w2 - x0 + x2);

        // Store new calculated vertices here.
        vec4 newPoint_10;
        vec4 newPoint_20;
        newPoint_10 = triangles[i].v1 + t_10*(triangles[i].v0-triangles[i].v1);
        newPoint_20 = triangles[i].v2 + t_20*(triangles[i].v0-triangles[i].v2);

        // Set outside vertex to one of the new vertices.
        triangles[i].v0 = newPoint_10;

        // Create a new triangle with the new vertices and one of the inside vertex.
        // Set the normals to be the same.
        Triangle extraTriangle(newPoint_10, newPoint_20, triangles[i].v2, triangles[i].color);
        extraTriangle.normal = triangles[i].normal;

        clippedTriangles.push_back(triangles[i]);
        clippedTriangles.push_back(extraTriangle);
      }
    }
    break;

    // Clip at right side x-plane
    case 2:
    for (int i = 0; i < triangles.size(); ++i) {
      dot[0] = triangles[i].v0.w * SCREEN_WIDTH/2;
      dot[1] = triangles[i].v1.w * SCREEN_WIDTH/2;
      dot[2] = triangles[i].v2.w * SCREEN_WIDTH/2;

      // CASE ALL VERTICES ARE IN THE PLANE
      if (triangles[i].v0.x < dot[0] && triangles[i].v1.x < dot[1] && triangles[i].v2.x < dot[2]) {
        clippedTriangles.push_back(triangles[i]);
      }

      // CASE V0 IS IN
      else if (triangles[i].v0.x < dot[0] && triangles[i].v1.x >= dot[1] && triangles[i].v2.x >= dot[2]) {
        // Get all x and w values.
        float x0 = triangles[i].v0.x;
        float x1 = triangles[i].v1.x;
        float x2 = triangles[i].v2.x;

        float w0 = triangles[i].v0.w;
        float w1 = triangles[i].v1.w;
        float w2 = triangles[i].v2.w;

        // Calculate the line intersection with the plane to change vertex positions within the plane.
        float t_01 = (x0 - SCREEN_WIDTH/2*w0)/(SCREEN_WIDTH/2*w1 - SCREEN_WIDTH/2*w0 - x1 + x0);
        float t_02 = (x0 - SCREEN_WIDTH/2*w0)/(SCREEN_WIDTH/2*w2 - SCREEN_WIDTH/2*w0 - x2 + x0);

        // Set new vertex positions.
        triangles[i].v1 = triangles[i].v0 + t_01*(triangles[i].v1-triangles[i].v0);
        triangles[i].v2 = triangles[i].v0 + t_02*(triangles[i].v2-triangles[i].v0);

        // Add new triangle to vector.
        clippedTriangles.push_back(triangles[i]);
      }

      // CASE V1 IS IN
      else if (triangles[i].v0.x >= dot[0] && triangles[i].v1.x < dot[1] && triangles[i].v2.x >= dot[2]) {
        // Get all x and w values.
        float x0 = triangles[i].v0.x;
        float x1 = triangles[i].v1.x;
        float x2 = triangles[i].v2.x;

        float w0 = triangles[i].v0.w;
        float w1 = triangles[i].v1.w;
        float w2 = triangles[i].v2.w;

        // Calculate the line intersection with the plane to change vertex positions within the plane.
        float t_10 = (x1 - SCREEN_WIDTH/2*w1)/(SCREEN_WIDTH/2*w0 - SCREEN_WIDTH/2*w1 - x0 + x1);
        float t_12 = (x1 - SCREEN_WIDTH/2*w1)/(SCREEN_WIDTH/2*w2 - SCREEN_WIDTH/2*w1 - x2 + x1);

        // Set new vertex positions.
        triangles[i].v0 = triangles[i].v1 + t_10*(triangles[i].v0-triangles[i].v1);
        triangles[i].v2 = triangles[i].v1 + t_12*(triangles[i].v2-triangles[i].v1);

        // Add new triangle to vector.
        clippedTriangles.push_back(triangles[i]);
      }

      // CASE V2 IS IN
      else if (triangles[i].v0.x >= dot[0] && triangles[i].v1.x >= dot[1] && triangles[i].v2.x < dot[2]) {
        // Get all x and w values.
        float x0 = triangles[i].v0.x;
        float x1 = triangles[i].v1.x;
        float x2 = triangles[i].v2.x;

        float w0 = triangles[i].v0.w;
        float w1 = triangles[i].v1.w;
        float w2 = triangles[i].v2.w;

        // Calculate the line intersection with the plane to change vertex positions within the plane.
        float t_21 = (x2 - SCREEN_WIDTH/2*w2)/(SCREEN_WIDTH/2*w1 - SCREEN_WIDTH/2*w2 - x1 + x2);
        float t_20 = (x2 - SCREEN_WIDTH/2*w2)/(SCREEN_WIDTH/2*w0 - SCREEN_WIDTH/2*w2 - x0 + x2);

        // Set new vertex positions.
        triangles[i].v1 = triangles[i].v2 + t_21*(triangles[i].v1-triangles[i].v2);
        triangles[i].v0 = triangles[i].v2 + t_20*(triangles[i].v0-triangles[i].v2);

        // Add new triangle to vector.
        clippedTriangles.push_back(triangles[i]);
      }

              /* WHEN 2 VERTICES ARE INSIDE, SPLIT INTO 2 TRIANGLES */

      // CASE V0 AND V1 ARE IN
      else if (triangles[i].v0.x < dot[0] && triangles[i].v1.x < dot[1] && triangles[i].v2.x >= dot[2]) {
        // Get all x and w values.
        float x0 = triangles[i].v0.x;
        float x1 = triangles[i].v1.x;
        float x2 = triangles[i].v2.x;

        float w0 = triangles[i].v0.w;
        float w1 = triangles[i].v1.w;
        float w2 = triangles[i].v2.w;

        // Calculate the line intersection with the plane to change vertex positions within the plane.
        float t_12 = (x1 - SCREEN_WIDTH/2*w1)/(SCREEN_WIDTH/2*w2 - SCREEN_WIDTH/2*w1 - x2 + x1);
        float t_02 = (x0 - SCREEN_WIDTH/2*w0)/(SCREEN_WIDTH/2*w2 - SCREEN_WIDTH/2*w0 - x2 + x0);

        // Store new calculated vertices here.
        vec4 newPoint_12;
        vec4 newPoint_02;
        newPoint_12 = triangles[i].v1 + t_12*(triangles[i].v2-triangles[i].v1);
        newPoint_02 = triangles[i].v0 + t_02*(triangles[i].v2-triangles[i].v0);

        // Set outside vertex to one of the new vertices.
        triangles[i].v2 = newPoint_02;

        // Create a new triangle with the new vertices and one of the inside vertex.
        // Set the normals to be the same.
        Triangle extraTriangle(newPoint_02, newPoint_12, triangles[i].v1, triangles[i].color);
        extraTriangle.normal = triangles[i].normal;

        // Add the new triangles to vector.
        clippedTriangles.push_back(triangles[i]);
        clippedTriangles.push_back(extraTriangle);
      }

      // CASE V0 AND V2 ARE IN
      else if (triangles[i].v0.x < dot[0] && triangles[i].v1.x >= dot[1] && triangles[i].v2.x < dot[2]) {
        // Get all x and w values.
        float x0 = triangles[i].v0.x;
        float x1 = triangles[i].v1.x;
        float x2 = triangles[i].v2.x;

        float w0 = triangles[i].v0.w;
        float w1 = triangles[i].v1.w;
        float w2 = triangles[i].v2.w;

        // Calculate the line intersection with the plane to change vertex positions within the plane.
        float t_01 = (x0 - SCREEN_WIDTH/2*w0)/(SCREEN_WIDTH/2*w1 - SCREEN_WIDTH/2*w0 - x1 + x0);
        float t_21 = (x2 - SCREEN_WIDTH/2*w2)/(SCREEN_WIDTH/2*w1 - SCREEN_WIDTH/2*w2 - x1 + x2);

        // Store new calculated vertices here.
        vec4 newPoint_01;
        vec4 newPoint_21;
        newPoint_01 = triangles[i].v0 + t_01*(triangles[i].v1-triangles[i].v0);
        newPoint_21 = triangles[i].v2 + t_21*(triangles[i].v1-triangles[i].v2);

        // Set outside vertex to one of the new vertices.
        triangles[i].v1 = newPoint_01;

        // Create a new triangle with the new vertices and one of the inside vertex.
        // Set the normals to be the same.
        Triangle extraTriangle(newPoint_01, newPoint_21, triangles[i].v2, triangles[i].color);
        extraTriangle.normal = triangles[i].normal;

        // Add the new triangles to vector.
        clippedTriangles.push_back(triangles[i]);
        clippedTriangles.push_back(extraTriangle);
      }

      // CASE V1 AND V2 ARE IN
      else if (triangles[i].v0.x >= dot[0] && triangles[i].v1.x < dot[1] && triangles[i].v2.x < dot[2]) {
        // Get all x and w values.
        float x0 = triangles[i].v0.x;
        float x1 = triangles[i].v1.x;
        float x2 = triangles[i].v2.x;

        float w0 = triangles[i].v0.w;
        float w1 = triangles[i].v1.w;
        float w2 = triangles[i].v2.w;

        // Calculate the line intersection with the plane to change vertex positions within the plane.
        float t_10 = (x1 - SCREEN_WIDTH/2*w1)/(SCREEN_WIDTH/2*w0 - SCREEN_WIDTH/2*w1 - x0 + x1);
        float t_20 = (x2 - SCREEN_WIDTH/2*w2)/(SCREEN_WIDTH/2*w0 - SCREEN_WIDTH/2*w2 - x0 + x2);

        // Store new calculated vertices here.
        vec4 newPoint_10;
        vec4 newPoint_20;
        newPoint_10 = triangles[i].v1 + t_10*(triangles[i].v0-triangles[i].v1);
        newPoint_20 = triangles[i].v2 + t_20*(triangles[i].v0-triangles[i].v2);

        // Set outside vertex to one of the new vertices.
        triangles[i].v0 = newPoint_10;

        // Create a new triangle with the new vertices and one of the inside vertex.
        // Set the normals to be the same.
        Triangle extraTriangle(newPoint_10, newPoint_20, triangles[i].v2, triangles[i].color);
        extraTriangle.normal = triangles[i].normal;

        // Add the new triangles to vector.
        clippedTriangles.push_back(triangles[i]);
        clippedTriangles.push_back(extraTriangle);
      }
    }
    break;

    // Clip at bottom side y-plane
    case 3:
    for (int i = 0; i < triangles.size(); ++i) {
      dot[0] = triangles[i].v0.w * SCREEN_HEIGHT/2;
      dot[1] = triangles[i].v1.w * SCREEN_HEIGHT/2;
      dot[2] = triangles[i].v2.w * SCREEN_HEIGHT/2;

      // CASE ALL VERTICES ARE IN THE PLANE
      if (triangles[i].v0.y < dot[0] && triangles[i].v1.y < dot[1] && triangles[i].v2.y < dot[2]) {
        clippedTriangles.push_back(triangles[i]);
      }

      // CASE V0 IS IN
      else if (triangles[i].v0.y < dot[0] && triangles[i].v1.y >= dot[1] && triangles[i].v2.y >= dot[2]) {
        // Get all x and w values.
        float y0 = triangles[i].v0.y;
        float y1 = triangles[i].v1.y;
        float y2 = triangles[i].v2.y;

        float w0 = triangles[i].v0.w;
        float w1 = triangles[i].v1.w;
        float w2 = triangles[i].v2.w;

        // Calculate the line intersection with the plane to change vertex positions within the plane.
        float t_01 = (y0 - SCREEN_HEIGHT/2*w0)/(SCREEN_HEIGHT/2*w1 - SCREEN_HEIGHT/2*w0 - y1 + y0);
        float t_02 = (y0 - SCREEN_HEIGHT/2*w0)/(SCREEN_HEIGHT/2*w2 - SCREEN_HEIGHT/2*w0 - y2 + y0);

        // Set new vertex positions.
        triangles[i].v1 = triangles[i].v0 + t_01*(triangles[i].v1-triangles[i].v0);
        triangles[i].v2 = triangles[i].v0 + t_02*(triangles[i].v2-triangles[i].v0);

        // Add new triangle to vector.
        clippedTriangles.push_back(triangles[i]);
      }

      // CASE V1 IS IN
      else if (triangles[i].v0.y >= dot[0] && triangles[i].v1.y < dot[1] && triangles[i].v2.y >= dot[2]) {
        // Get all x and w values.
        float y0 = triangles[i].v0.y;
        float y1 = triangles[i].v1.y;
        float y2 = triangles[i].v2.y;

        float w0 = triangles[i].v0.w;
        float w1 = triangles[i].v1.w;
        float w2 = triangles[i].v2.w;

        // Calculate the line intersection with the plane to change vertex positions within the plane.
        float t_10 = (y1 - SCREEN_HEIGHT/2*w1)/(SCREEN_HEIGHT/2*w0 - SCREEN_HEIGHT/2*w1 - y0 + y1);
        float t_12 = (y1 - SCREEN_HEIGHT/2*w1)/(SCREEN_HEIGHT/2*w2 - SCREEN_HEIGHT/2*w1 - y2 + y1);

        // Set new vertex positions.
        triangles[i].v0 = triangles[i].v1 + t_10*(triangles[i].v0-triangles[i].v1);
        triangles[i].v2 = triangles[i].v1 + t_12*(triangles[i].v2-triangles[i].v1);

        // Add new triangle to vector.
        clippedTriangles.push_back(triangles[i]);
      }

      // CASE V2 IS IN
      else if (triangles[i].v0.y >= dot[0] && triangles[i].v1.y >= dot[1] && triangles[i].v2.y < dot[2]) {
        // Get all x and w values.
        float y0 = triangles[i].v0.y;
        float y1 = triangles[i].v1.y;
        float y2 = triangles[i].v2.y;

        float w0 = triangles[i].v0.w;
        float w1 = triangles[i].v1.w;
        float w2 = triangles[i].v2.w;

        // Calculate the line intersection with the plane to change vertex positions within the plane.
        float t_21 = (y2 - SCREEN_HEIGHT/2*w2)/(SCREEN_HEIGHT/2*w1 - SCREEN_HEIGHT/2*w2 - y1 + y2);
        float t_20 = (y2 - SCREEN_HEIGHT/2*w2)/(SCREEN_HEIGHT/2*w0 - SCREEN_HEIGHT/2*w2 - y0 + y2);

        triangles[i].v1 = triangles[i].v2 + t_21*(triangles[i].v1-triangles[i].v2);
        triangles[i].v0 = triangles[i].v2 + t_20*(triangles[i].v0-triangles[i].v2);

        // Add new triangle to vector.
        clippedTriangles.push_back(triangles[i]);
      }

              /* WHEN 2 VERTICES ARE INSIDE, SPLIT INTO 2 TRIANGLES */

      // CASE V0 AND V1 ARE IN
      else if (triangles[i].v0.y < dot[0] && triangles[i].v1.y < dot[1] && triangles[i].v2.y >= dot[2]) {
        // Get all x and w values.
        float y0 = triangles[i].v0.y;
        float y1 = triangles[i].v1.y;
        float y2 = triangles[i].v2.y;

        float w0 = triangles[i].v0.w;
        float w1 = triangles[i].v1.w;
        float w2 = triangles[i].v2.w;

        // Calculate the line intersection with the plane to change vertex positions within the plane.
        float t_12 = (y1 - SCREEN_HEIGHT/2*w1)/(SCREEN_HEIGHT/2*w2 - SCREEN_HEIGHT/2*w1 - y2 + y1);
        float t_02 = (y0 - SCREEN_HEIGHT/2*w0)/(SCREEN_HEIGHT/2*w2 - SCREEN_HEIGHT/2*w0 - y2 + y0);

        // Store new calculated vertices here.
        vec4 newPoint_12;
        vec4 newPoint_02;
        newPoint_12 = triangles[i].v1 + t_12*(triangles[i].v2-triangles[i].v1);
        newPoint_02 = triangles[i].v0 + t_02*(triangles[i].v2-triangles[i].v0);

        // Set outside vertex to one of the new vertices.
        triangles[i].v2 = newPoint_02;

        // Create a new triangle with the new vertices and one of the inside vertex.
        // Set the normals to be the same.
        Triangle extraTriangle(newPoint_02, newPoint_12, triangles[i].v1, triangles[i].color);
        extraTriangle.normal = triangles[i].normal;

        // Add the new triangles to vector.
        clippedTriangles.push_back(triangles[i]);
        clippedTriangles.push_back(extraTriangle);
      }

      // CASE V0 AND V2 ARE IN
      else if (triangles[i].v0.y < dot[0] && triangles[i].v1.y >= dot[1] && triangles[i].v2.y < dot[2]) {
        // Get all x and w values.
        float y0 = triangles[i].v0.y;
        float y1 = triangles[i].v1.y;
        float y2 = triangles[i].v2.y;

        float w0 = triangles[i].v0.w;
        float w1 = triangles[i].v1.w;
        float w2 = triangles[i].v2.w;

        // Calculate the line intersection with the plane to change vertex positions within the plane.
        float t_01 = (y0 - SCREEN_HEIGHT/2*w0)/(SCREEN_HEIGHT/2*w1 - SCREEN_HEIGHT/2*w0 - y1 + y0);
        float t_21 = (y2 - SCREEN_HEIGHT/2*w2)/(SCREEN_HEIGHT/2*w1 - SCREEN_HEIGHT/2*w2 - y1 + y2);

        // Store new calculated vertices here.
        vec4 newPoint_01;
        vec4 newPoint_21;
        newPoint_01 = triangles[i].v0 + t_01*(triangles[i].v1-triangles[i].v0);
        newPoint_21 = triangles[i].v2 + t_21*(triangles[i].v1-triangles[i].v2);

        // Set outside vertex to one of the new vertices.
        triangles[i].v1 = newPoint_01;

        // Create a new triangle with the new vertices and one of the inside vertex.
        // Set the normals to be the same.
        Triangle extraTriangle(newPoint_01, newPoint_21, triangles[i].v2, triangles[i].color);
        extraTriangle.normal = triangles[i].normal;

        // Add the new triangles to vector.
        clippedTriangles.push_back(triangles[i]);
        clippedTriangles.push_back(extraTriangle);
      }

      // CASE V1 AND V2 ARE IN
      else if (triangles[i].v0.y >= dot[0] && triangles[i].v1.y < dot[1] && triangles[i].v2.y < dot[2]) {
        // Get all x and w values.
        float y0 = triangles[i].v0.y;
        float y1 = triangles[i].v1.y;
        float y2 = triangles[i].v2.y;

        float w0 = triangles[i].v0.w;
        float w1 = triangles[i].v1.w;
        float w2 = triangles[i].v2.w;

        // Calculate the line intersection with the plane to change vertex positions within the plane.
        float t_10 = (y1 - SCREEN_HEIGHT/2*w1)/(SCREEN_HEIGHT/2*w0 - SCREEN_HEIGHT/2*w1 - y0 + y1);
        float t_20 = (y2 - SCREEN_HEIGHT/2*w2)/(SCREEN_HEIGHT/2*w0 - SCREEN_HEIGHT/2*w2 - y0 + y2);

        // Store new calculated vertices here.
        vec4 newPoint_10;
        vec4 newPoint_20;
        newPoint_10 = triangles[i].v1 + t_10*(triangles[i].v0-triangles[i].v1);
        newPoint_20 = triangles[i].v2 + t_20*(triangles[i].v0-triangles[i].v2);

        // Set outside vertex to one of the new vertices.
        triangles[i].v0 = newPoint_10;

        // Create a new triangle with the new vertices and one of the inside vertex.
        // Set the normals to be the same.
        Triangle extraTriangle(newPoint_10, newPoint_20, triangles[i].v2, triangles[i].color);
        extraTriangle.normal = triangles[i].normal;

        // Add the new triangles to vector.
        clippedTriangles.push_back(triangles[i]);
        clippedTriangles.push_back(extraTriangle);
      }
    }
    break;

    // Clip at top side y-plane
    case 4:
    for (int i = 0; i < triangles.size(); ++i) {
      dot[0] = triangles[i].v0.w * -SCREEN_HEIGHT/2;
      dot[1] = triangles[i].v1.w * -SCREEN_HEIGHT/2;
      dot[2] = triangles[i].v2.w * -SCREEN_HEIGHT/2;

      // CASE ALL VERTICES ARE IN THE PLANE
      if (triangles[i].v0.y > dot[0] && triangles[i].v1.y > dot[1] && triangles[i].v2.y > dot[2]) {
        clippedTriangles.push_back(triangles[i]);
      }

      // CASE V0 IS IN
      else if (triangles[i].v0.y > dot[0] && triangles[i].v1.y <= dot[1] && triangles[i].v2.y <= dot[2]) {
        // Get all x and w values.
        float y0 = triangles[i].v0.y;
        float y1 = triangles[i].v1.y;
        float y2 = triangles[i].v2.y;

        float w0 = triangles[i].v0.w;
        float w1 = triangles[i].v1.w;
        float w2 = triangles[i].v2.w;

        // Calculate the line intersection with the plane to change vertex positions within the plane.
        float t_01 = (y0 + SCREEN_HEIGHT/2*w0)/(-SCREEN_HEIGHT/2*w1 + SCREEN_HEIGHT/2*w0 - y1 + y0);
        float t_02 = (y0 + SCREEN_HEIGHT/2*w0)/(-SCREEN_HEIGHT/2*w2 + SCREEN_HEIGHT/2*w0 - y2 + y0);

        // Set new vertex positions.
        triangles[i].v1 = triangles[i].v0 + t_01*(triangles[i].v1-triangles[i].v0);
        triangles[i].v2 = triangles[i].v0 + t_02*(triangles[i].v2-triangles[i].v0);

        // Add new triangle to vector.
        clippedTriangles.push_back(triangles[i]);
      }

      // CASE V1 IS IN
      else if (triangles[i].v0.y <= dot[0] && triangles[i].v1.y > dot[1] && triangles[i].v2.y <= dot[2]) {
        // Get all x and w values.
        float y0 = triangles[i].v0.y;
        float y1 = triangles[i].v1.y;
        float y2 = triangles[i].v2.y;

        float w0 = triangles[i].v0.w;
        float w1 = triangles[i].v1.w;
        float w2 = triangles[i].v2.w;

        // Calculate the line intersection with the plane to change vertex positions within the plane.
        float t_10 = (y1 + SCREEN_HEIGHT/2*w1)/(-SCREEN_HEIGHT/2*w0 + SCREEN_HEIGHT/2*w1 - y0 + y1);
        float t_12 = (y1 + SCREEN_HEIGHT/2*w1)/(-SCREEN_HEIGHT/2*w2 + SCREEN_HEIGHT/2*w1 - y2 + y1);

        // Set new vertex positions.
        triangles[i].v0 = triangles[i].v1 + t_10*(triangles[i].v0-triangles[i].v1);
        triangles[i].v2 = triangles[i].v1 + t_12*(triangles[i].v2-triangles[i].v1);

        // Add new triangle to vector.
        clippedTriangles.push_back(triangles[i]);
      }

      // CASE V2 IS IN
      else if (triangles[i].v0.y <= dot[0] && triangles[i].v1.y <= dot[1] && triangles[i].v2.y > dot[2]) {
        // Get all x and w values.
        float y0 = triangles[i].v0.y;
        float y1 = triangles[i].v1.y;
        float y2 = triangles[i].v2.y;

        float w0 = triangles[i].v0.w;
        float w1 = triangles[i].v1.w;
        float w2 = triangles[i].v2.w;

        // Calculate the line intersection with the plane to change vertex positions within the plane.
        float t_21 = (y2 + SCREEN_HEIGHT/2*w2)/(-SCREEN_HEIGHT/2*w1 + SCREEN_HEIGHT/2*w2 - y1 + y2);
        float t_20 = (y2 + SCREEN_HEIGHT/2*w2)/(-SCREEN_HEIGHT/2*w0 + SCREEN_HEIGHT/2*w2 - y0 + y2);

        // Set new vertex positions.
        triangles[i].v1 = triangles[i].v2 + t_21*(triangles[i].v1-triangles[i].v2);
        triangles[i].v0 = triangles[i].v2 + t_20*(triangles[i].v0-triangles[i].v2);

        // Add new triangle to vector.
        clippedTriangles.push_back(triangles[i]);
      }
              /* WHEN 2 VERTICES ARE INSIDE, SPLIT INTO 2 TRIANGLES */

      // CASE V0 AND V1 ARE IN
      else if (triangles[i].v0.y > dot[0] && triangles[i].v1.y > dot[1] && triangles[i].v2.y <= dot[2]) {
        // Get all x and w values.
        float y0 = triangles[i].v0.y;
        float y1 = triangles[i].v1.y;
        float y2 = triangles[i].v2.y;

        float w0 = triangles[i].v0.w;
        float w1 = triangles[i].v1.w;
        float w2 = triangles[i].v2.w;

        // Calculate the line intersection with the plane to change vertex positions within the plane.
        float t_12 = (y1 + SCREEN_HEIGHT/2*w1)/(-SCREEN_HEIGHT/2*w2 + SCREEN_HEIGHT/2*w1 - y2 + y1);
        float t_02 = (y0 + SCREEN_HEIGHT/2*w0)/(-SCREEN_HEIGHT/2*w2 + SCREEN_HEIGHT/2*w0 - y2 + y0);

        // Store new calculated vertices here.
        vec4 newPoint_12;
        vec4 newPoint_02;
        newPoint_12 = triangles[i].v1 + t_12*(triangles[i].v2-triangles[i].v1);
        newPoint_02 = triangles[i].v0 + t_02*(triangles[i].v2-triangles[i].v0);

        // Set outside vertex to one of the new vertices.
        triangles[i].v2 = newPoint_02;

        // Create a new triangle with the new vertices and one of the inside vertex.
        // Set the normals to be the same.
        Triangle extraTriangle(newPoint_02, newPoint_12, triangles[i].v1, triangles[i].color);
        extraTriangle.normal = triangles[i].normal;

        // Add the new triangles to vector.
        clippedTriangles.push_back(triangles[i]);
        clippedTriangles.push_back(extraTriangle);
      }

      // CASE V0 AND V2 ARE IN
      else if (triangles[i].v0.y > dot[0] && triangles[i].v1.y <= dot[1] && triangles[i].v2.y > dot[2]) {
        // Get all x and w values.
        float y0 = triangles[i].v0.y;
        float y1 = triangles[i].v1.y;
        float y2 = triangles[i].v2.y;

        float w0 = triangles[i].v0.w;
        float w1 = triangles[i].v1.w;
        float w2 = triangles[i].v2.w;

        // Calculate the line intersection with the plane to change vertex positions within the plane.
        float t_01 = (y0 + SCREEN_HEIGHT/2*w0)/(-SCREEN_HEIGHT/2*w1 + SCREEN_HEIGHT/2*w0 - y1 + y0);
        float t_21 = (y2 + SCREEN_HEIGHT/2*w2)/(-SCREEN_HEIGHT/2*w1 + SCREEN_HEIGHT/2*w2 - y1 + y2);

        // Store new calculated vertices here.
        vec4 newPoint_01;
        vec4 newPoint_21;
        newPoint_01 = triangles[i].v0 + t_01*(triangles[i].v1-triangles[i].v0);
        newPoint_21 = triangles[i].v2 + t_21*(triangles[i].v1-triangles[i].v2);

        // Set outside vertex to one of the new vertices.
        triangles[i].v1 = newPoint_01;

        // Create a new triangle with the new vertices and one of the inside vertex.
        // Set the normals to be the same.
        Triangle extraTriangle(newPoint_01, newPoint_21, triangles[i].v2, triangles[i].color);
        extraTriangle.normal = triangles[i].normal;

        // Add the new triangles to vector.
        clippedTriangles.push_back(triangles[i]);
        clippedTriangles.push_back(extraTriangle);
      }

      // CASE V1 AND V2 ARE IN
      else if (triangles[i].v0.y <= dot[0] && triangles[i].v1.y > dot[1] && triangles[i].v2.y > dot[2]) {
        // Get all x and w values.
        float y0 = triangles[i].v0.y;
        float y1 = triangles[i].v1.y;
        float y2 = triangles[i].v2.y;

        float w0 = triangles[i].v0.w;
        float w1 = triangles[i].v1.w;
        float w2 = triangles[i].v2.w;

        // Calculate the line intersection with the plane to change vertex positions within the plane.
        float t_10 = (y1 + SCREEN_HEIGHT/2*w1)/(-SCREEN_HEIGHT/2*w0 + SCREEN_HEIGHT/2*w1 - y0 + y1);
        float t_20 = (y2 + SCREEN_HEIGHT/2*w2)/(-SCREEN_HEIGHT/2*w0 + SCREEN_HEIGHT/2*w2 - y0 + y2);

        // Store new calculated vertices here.
        vec4 newPoint_10;
        vec4 newPoint_20;
        newPoint_10 = triangles[i].v1 + t_10*(triangles[i].v0-triangles[i].v1);
        newPoint_20 = triangles[i].v2 + t_20*(triangles[i].v0-triangles[i].v2);

        // Set outside vertex to one of the new vertices.
        triangles[i].v0 = newPoint_10;

        // Create a new triangle with the new vertices and one of the inside vertex.
        // Set the normals to be the same.
        Triangle extraTriangle(newPoint_10, newPoint_20, triangles[i].v2, triangles[i].color);
        extraTriangle.normal = triangles[i].normal;

        // Add the new triangles to vector.
        clippedTriangles.push_back(triangles[i]);
        clippedTriangles.push_back(extraTriangle);
      }
    }
    break;

    // Clip at z plane at current camera position (anything behind the camera)
    case 5:
    for (int i = 0; i < triangles.size(); ++i) {

      // CASE ALL VERTICES ARE IN THE PLANE
      if (triangles[i].v0.z > 0.01f && triangles[i].v1.z > 0.01f && triangles[i].v2.z > 0.01f) {
        clippedTriangles.push_back(triangles[i]);
      }

      // // CASE V0 IS IN
      // else if (triangles[i].v0.z > 0.01f && triangles[i].v1.z <= 0.01f && triangles[i].v2.z <= 0.01f) {
      //   // Get w values.
      //   float z0 = triangles[i].v0.z;
      //   float z1 = triangles[i].v1.z;
      //   float z2 = triangles[i].v2.z;
      //
      //   // Calculate the line intersection with the plane to change vertex positions within the plane.
      //   float t_01 = (0.01f-z0)/(z1 - z0);
      //   printf("%f\n", t_01);
      //   float t_02 = (0.01f-z0)/(z2 - z0);
      //   printf("%f\n", t_02);
      //
      //   // Set new vertex positions.
      //   printf("%f, %f, %f\n", triangles[i].v1.z, triangles[i].v0.z, t_01);
      //   printf("\n" );
      //   triangles[i].v1 = triangles[i].v0 + t_01*(triangles[i].v1-triangles[i].v0);
      //   triangles[i].v2 = triangles[i].v0 + t_02*(triangles[i].v2-triangles[i].v0);
      //   printf("%f, %f, %f, %f\n", triangles[i].v0.x, triangles[i].v0.y, triangles[i].v0.z, triangles[i].v0.w);
      //   printf("%f, %f, %f, %f\n", triangles[i].v1.x, triangles[i].v1.y, triangles[i].v1.z, triangles[i].v1.w);
      //   printf("%f, %f, %f, %f\n", triangles[i].v2.x, triangles[i].v2.y, triangles[i].v2.z, triangles[i].v2.w);
      //   printf("\n" );
      //
      //   // Add new triangle to vector.
      //   clippedTriangles.push_back(triangles[i]);
      // }

      // // CASE V1 IS IN
      // else if (triangles[i].v0.w <= 0.01f && triangles[i].v1.w > 0.01f && triangles[i].v2.w <= 0.01f) {
      //   // Get w values.
      //   float w0 = triangles[i].v0.w;
      //   float w1 = triangles[i].v1.w;
      //   float w2 = triangles[i].v2.w;
      //
      //   // Calculate the line intersection with the plane to change vertex positions within the plane.
      //   float t_10 = w1/(w1 - w0);
      //   float t_12 = w1/(w1 - w2);
      //
      //   // Set new vertex positions.
      //   triangles[i].v0 = triangles[i].v1 + t_10*(triangles[i].v0-triangles[i].v1);
      //   triangles[i].v2 = triangles[i].v1 + t_12*(triangles[i].v2-triangles[i].v1);
      //
      //   // Add new triangle to vector.
      //   clippedTriangles.push_back(triangles[i]);
      // }

      // // CASE V2 IS IN
      // else if (triangles[i].v0.w <= 0.01f && triangles[i].v1.w <= 0.01f && triangles[i].v2.w > 0.01f) {
      //   // Get w values.
      //   float w0 = triangles[i].v0.w;
      //   float w1 = triangles[i].v1.w;
      //   float w2 = triangles[i].v2.w;
      //
      //   // Calculate the line intersection with the plane to change vertex positions within the plane.
      //   float t_21 = w2/(w2 - w1);
      //   float t_20 = w2/(w2 - w0);
      //
      //   // Set new vertex positions.
      //   triangles[i].v1 = triangles[i].v2 + t_21*(triangles[i].v1-triangles[i].v2);
      //   triangles[i].v0 = triangles[i].v2 + t_20*(triangles[i].v0-triangles[i].v2);
      //
      //   // Add new triangle to vector.
      //   clippedTriangles.push_back(triangles[i]);
      // }
      //         /* WHEN 2 VERTICES ARE INSIDE, SPLIT INTO 2 TRIANGLES */
      //
      // // CASE V0 AND V1 ARE IN
      // else if (triangles[i].v0.w > 0.01f && triangles[i].v1.w > 0.01f && triangles[i].v2.w <= 0.01f) {
      //   // Get w values.
      //   float w0 = triangles[i].v0.w;
      //   float w1 = triangles[i].v1.w;
      //   float w2 = triangles[i].v2.w;
      //
      //   // Calculate the line intersection with the plane to change vertex positions within the plane.
      //   float t_12 = w1/(w1 - w2);
      //   float t_02 = w0/(w0 - w2);
      //
      //   // Store new calculated vertices here.
      //   vec4 newPoint_12;
      //   vec4 newPoint_02;
      //   newPoint_12 = triangles[i].v1 + t_12*(triangles[i].v2-triangles[i].v1);
      //   newPoint_02 = triangles[i].v0 + t_02*(triangles[i].v2-triangles[i].v0);
      //
      //   // Set outside vertex to one of the new vertices.
      //   triangles[i].v2 = newPoint_02;
      //
      //   // Create a new triangle with the new vertices and one of the inside vertex.
      //   // Set the normals to be the same.
      //   Triangle extraTriangle(newPoint_02, newPoint_12, triangles[i].v1, triangles[i].color);
      //   extraTriangle.normal = triangles[i].normal;
      //
      //   // Add the new triangles to vector.
      //   clippedTriangles.push_back(triangles[i]);
      //   clippedTriangles.push_back(extraTriangle);
      // }
      //
      // // CASE V0 AND V2 ARE IN
      // else if (triangles[i].v0.w > 0.01f && triangles[i].v1.w <= 0.01f && triangles[i].v2.x > 0.01f) {
      //   printf("YES\n" );
      //   // Get w values.
      //   float w0 = triangles[i].v0.w;
      //   float w1 = triangles[i].v1.w;
      //   float w2 = triangles[i].v2.w;
      //
      //   // Calculate the line intersection with the plane to change vertex positions within the plane.
      //   float t_01 = w0/(w0 - w1);
      //   float t_21 = w2/(w2 - w1);
      //
      //   // Store new calculated vertices here.
      //   vec4 newPoint_01;
      //   vec4 newPoint_21;
      //   newPoint_01 = triangles[i].v0 + t_01*(triangles[i].v1-triangles[i].v0);
      //   newPoint_21 = triangles[i].v2 + t_21*(triangles[i].v1-triangles[i].v2);
      //
      //   // Set outside vertex to one of the new vertices.
      //   triangles[i].v1 = newPoint_01;
      //
      //   // Create a new triangle with the new vertices and one of the inside vertex.
      //   // Set the normals to be the same.
      //   Triangle extraTriangle(newPoint_01, newPoint_21, triangles[i].v2, triangles[i].color);
      //   extraTriangle.normal = triangles[i].normal;
      //
      //   // Add the new triangles to vector.
      //   clippedTriangles.push_back(triangles[i]);
      //   clippedTriangles.push_back(extraTriangle);
      // }
      //
      // // CASE V1 AND V2 ARE IN
      // else if (triangles[i].v0.w <= 0.01f && triangles[i].v1.w > 0.01f && triangles[i].v2.w > 0.01f) {
      //   // Get w values.
      //   float w0 = triangles[i].v0.w;
      //   float w1 = triangles[i].v1.w;
      //   float w2 = triangles[i].v2.w;
      //
      //   // Calculate the line intersection with the plane to change vertex positions within the plane.
      //   float t_10 = w1/(w1 - w0);
      //   float t_20 = w2/(w2 - w0);
      //
      //   // Store new calculated vertices here.
      //   vec4 newPoint_10;
      //   vec4 newPoint_20;
      //   newPoint_10 = triangles[i].v1 + t_10*(triangles[i].v0-triangles[i].v1);
      //   newPoint_20 = triangles[i].v2 + t_20*(triangles[i].v0-triangles[i].v2);
      //
      //   // Set outside vertex to one of the new vertices.
      //   triangles[i].v0 = newPoint_10;
      //
      //   // Create a new triangle with the new vertices and one of the inside vertex.
      //   // Set the normals to be the same.
      //   Triangle extraTriangle(newPoint_10, newPoint_20, triangles[i].v2, triangles[i].color);
      //   extraTriangle.normal = triangles[i].normal;
      //
      //   // Add the new triangles to vector.
      //   clippedTriangles.push_back(triangles[i]);
      //   clippedTriangles.push_back(extraTriangle);
      // }
    }
    break;
    // Clip at far z plane (now set to z to 10)
    case 6:
    for (int i = 0; i < triangles.size(); ++i) {
      float wlimit = 5.0f/focalLength;

      // CASE ALL VERTICES ARE IN THE PLANE
      if (triangles[i].v0.w <= wlimit && triangles[i].v1.w <= wlimit && triangles[i].v2.w <= wlimit) {
        clippedTriangles.push_back(triangles[i]);
      }

      // CASE V0 IS IN
      else if (triangles[i].v0.w <= wlimit && triangles[i].v1.w > wlimit && triangles[i].v2.w > wlimit) {
        // Get w values.
        float w0 = triangles[i].v0.w;
        float w1 = triangles[i].v1.w;
        float w2 = triangles[i].v2.w;

        // Calculate the line intersection with the plane to change vertex positions within the plane.
        float t_01 = (wlimit-w0)/(w1 - w0);
        float t_02 = (wlimit-w0)/(w2 - w0);

        // Set new vertex positions.
        triangles[i].v1 = triangles[i].v0 + t_01*(triangles[i].v1-triangles[i].v0);
        triangles[i].v2 = triangles[i].v0 + t_02*(triangles[i].v2-triangles[i].v0);

        // Add new triangle to vector.
        clippedTriangles.push_back(triangles[i]);
      }

      // CASE V1 IS IN
      else if (triangles[i].v0.w > wlimit && triangles[i].v1.w <= wlimit && triangles[i].v2.w > wlimit) {
        // Get w values.
        float w0 = triangles[i].v0.w;
        float w1 = triangles[i].v1.w;
        float w2 = triangles[i].v2.w;

        // Calculate the line intersection with the plane to change vertex positions within the plane.
        float t_10 = (wlimit-w1)/(w0 - w1);
        float t_12 = (wlimit-w1)/(w2 - w1);

        // Set new vertex positions.
        triangles[i].v0 = triangles[i].v1 + t_10*(triangles[i].v0-triangles[i].v1);
        triangles[i].v2 = triangles[i].v1 + t_12*(triangles[i].v2-triangles[i].v1);

        // Add new triangle to vector.
        clippedTriangles.push_back(triangles[i]);
      }

      // CASE V2 IS IN
      else if (triangles[i].v0.w > wlimit && triangles[i].v1.w > wlimit && triangles[i].v2.w <= wlimit) {
        // Get w values.
        float w0 = triangles[i].v0.w;
        float w1 = triangles[i].v1.w;
        float w2 = triangles[i].v2.w;

        // Calculate the line intersection with the plane to change vertex positions within the plane.
        float t_21 = (wlimit-w2)/(w1 - w2);
        float t_20 = (wlimit-w2)/(w0 - w2);

        // Set new vertex positions.
        triangles[i].v1 = triangles[i].v2 + t_21*(triangles[i].v1-triangles[i].v2);
        triangles[i].v0 = triangles[i].v2 + t_20*(triangles[i].v0-triangles[i].v2);

        // Add new triangle to vector.
        clippedTriangles.push_back(triangles[i]);
      }
              /* WHEN 2 VERTICES ARE INSIDE, SPLIT INTO 2 TRIANGLES */

      // CASE V0 AND V1 ARE IN
      else if (triangles[i].v0.w <= wlimit && triangles[i].v1.w <= wlimit && triangles[i].v2.w > wlimit) {
        // Get w values.
        float w0 = triangles[i].v0.w;
        float w1 = triangles[i].v1.w;
        float w2 = triangles[i].v2.w;

        // Calculate the line intersection with the plane to change vertex positions within the plane.
        float t_12 = (wlimit-w1)/(w2 - w1);
        float t_02 = (wlimit-w0)/(w2 - w0);

        // Store new calculated vertices here.
        vec4 newPoint_12;
        vec4 newPoint_02;
        newPoint_12 = triangles[i].v1 + t_12*(triangles[i].v2-triangles[i].v1);
        newPoint_02 = triangles[i].v0 + t_02*(triangles[i].v2-triangles[i].v0);

        // Set outside vertex to one of the new vertices.
        triangles[i].v2 = newPoint_02;

        // Create a new triangle with the new vertices and one of the inside vertex.
        // Set the normals to be the same.
        Triangle extraTriangle(newPoint_02, newPoint_12, triangles[i].v1, triangles[i].color);
        extraTriangle.normal = triangles[i].normal;

        // Add the new triangles to vector.
        clippedTriangles.push_back(triangles[i]);
        clippedTriangles.push_back(extraTriangle);
      }

      // CASE V0 AND V2 ARE IN
      else if (triangles[i].v0.w <= wlimit && triangles[i].v1.w > wlimit && triangles[i].v2.x <= wlimit) {
        // Get w values.
        float w0 = triangles[i].v0.w;
        float w1 = triangles[i].v1.w;
        float w2 = triangles[i].v2.w;

        // Calculate the line intersection with the plane to change vertex positions within the plane.
        float t_01 = (wlimit-w0)/(w1 - w0);
        float t_21 = (wlimit-w2)/(w1 - w0);

        // Store new calculated vertices here.
        vec4 newPoint_01;
        vec4 newPoint_21;
        newPoint_01 = triangles[i].v0 + t_01*(triangles[i].v1-triangles[i].v0);
        newPoint_21 = triangles[i].v2 + t_21*(triangles[i].v1-triangles[i].v2);

        // Set outside vertex to one of the new vertices.
        triangles[i].v1 = newPoint_01;

        // Create a new triangle with the new vertices and one of the inside vertex.
        // Set the normals to be the same.
        Triangle extraTriangle(newPoint_01, newPoint_21, triangles[i].v2, triangles[i].color);
        extraTriangle.normal = triangles[i].normal;

        // Add the new triangles to vector.
        clippedTriangles.push_back(triangles[i]);
        clippedTriangles.push_back(extraTriangle);
      }

      // CASE V1 AND V2 ARE IN
      else if (triangles[i].v0.w > wlimit && triangles[i].v1.w <= wlimit && triangles[i].v2.w <= wlimit) {
        // Get w values.
        float w0 = triangles[i].v0.w;
        float w1 = triangles[i].v1.w;
        float w2 = triangles[i].v2.w;

        // Calculate the line intersection with the plane to change vertex positions within the plane.
        float t_10 = (wlimit-w1)/(w0 - w1);
        float t_20 = (wlimit-w2)/(w0 - w2);

        // Store new calculated vertices here.
        vec4 newPoint_10;
        vec4 newPoint_20;
        newPoint_10 = triangles[i].v1 + t_10*(triangles[i].v0-triangles[i].v1);
        newPoint_20 = triangles[i].v2 + t_20*(triangles[i].v0-triangles[i].v2);

        // Set outside vertex to one of the new vertices.
        triangles[i].v0 = newPoint_10;

        // Create a new triangle with the new vertices and one of the inside vertex.
        // Set the normals to be the same.
        Triangle extraTriangle(newPoint_10, newPoint_20, triangles[i].v2, triangles[i].color);
        extraTriangle.normal = triangles[i].normal;

        // Add the new triangles to vector.
        clippedTriangles.push_back(triangles[i]);
        clippedTriangles.push_back(extraTriangle);
      }
    }
    break;
  }
  return clippedTriangles;
}

/* THIS FUNCTION ADDS NEW POLYGONS TO SIMULATE SHADOWS */
vector<Triangle> createShadowVolume( vector<Triangle>& triangles ) {

  // Vector to store all triangles with their shadow volumes.
  vector<Triangle> trianglesWithShadows;

  // Store new vertices here.
  vec4 n0, n1, n2;

  // Go through the vector of triangles.
  for (int i = 0; i < triangles.size(); ++i) {
    // Add original triangle.
    trianglesWithShadows.push_back(triangles[i]);

    // Original vertices.
    vec4 v0 = triangles[i].v0;
    vec4 v1 = triangles[i].v1;
    vec4 v2 = triangles[i].v2;

    // Calculate position of new vertices. (Light -> triangle vertex) * arbitrary number.
    n0 = ((triangles[i].v0 - (lightPos-cameraPos)) * 100.0f);
    n1 = ((triangles[i].v1 - (lightPos-cameraPos)) * 100.0f);
    n2 = ((triangles[i].v2 - (lightPos-cameraPos)) * 100.0f);

    // Triangle shadowVolume0(v0, v1, v2, vec3(1.0f,1.0f,1.0f));
    // shadowVolume0.normal = triangles[i].normal;
    // Triangle shadowVolume1(n0, n1, n2, vec3(1.0f,1.0f,1.0f));
    // shadowVolume1.normal = triangles[i].normal;

    // Connect the vertices together to form triangles.
    Triangle shadowVolume2(v0, n0, v1, vec3(-1.0f,-1.0f,-1.0f));
    Triangle shadowVolume3(n0, v1, n1, vec3(-1.0f,-1.0f,-1.0f));
    Triangle shadowVolume4(v1, n1, v2, vec3(-1.0f,-1.0f,-1.0f));
    Triangle shadowVolume5(n1, v2, n2, vec3(-1.0f,-1.0f,-1.0f));
    Triangle shadowVolume6(v2, n2, v0, vec3(-1.0f,-1.0f,-1.0f));
    Triangle shadowVolume7(n2, v0, n0, vec3(-1.0f,-1.0f,-1.0f));

    // Compute normals of all triangles.
    shadowVolume2.ComputeNormal();
    shadowVolume3.ComputeNormal();
    shadowVolume4.ComputeNormal();
    shadowVolume5.ComputeNormal();
    shadowVolume6.ComputeNormal();
    shadowVolume7.ComputeNormal();

    // Add shadow triangles to the vector.
    trianglesWithShadows.push_back(shadowVolume2);
    trianglesWithShadows.push_back(shadowVolume3);
    trianglesWithShadows.push_back(shadowVolume4);
    trianglesWithShadows.push_back(shadowVolume5);
    trianglesWithShadows.push_back(shadowVolume6);
    trianglesWithShadows.push_back(shadowVolume7);
  }

  return trianglesWithShadows;
}

float surroundingShadowSum(int y, int x) {
  float val = 0;
  val = shadowBuffer[y][x] + shadowBuffer[y-1][x] + shadowBuffer[y-1][x-1] + shadowBuffer[y-1][x+1] +
        shadowBuffer[y+1][x-1] + shadowBuffer[y+1][x] + shadowBuffer[y+1][x-1] + shadowBuffer[y][x-1] +
        shadowBuffer[y][x+1];
  val /= 9.0f;

  return val;
}
