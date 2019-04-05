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

#define SCREEN_WIDTH 320
#define SCREEN_HEIGHT 256
#define FULLSCREEN_MODE false


/* ----------------------------------------------------------------------------*/
/* VARIABLES                                                                   */
/* ----------------------------------------------------------------------------*/

float focalLength = 256;
vec4 cameraPos( 0, 0, -3.001, 1 );
// float yaw = 0.0;
// mat4 R(1.0f);
vec3 currentColor;
float depthBuffer[SCREEN_HEIGHT][SCREEN_WIDTH]; // Stores depth information of screen pixels. This
                                                // is stored as 1/z.

// Initialise light variables
vec4 lightPos(0,-0.5,-0.7, 1);
vec3 lightPower = 11.0f*vec3( 1, 1, 1 );
vec3 indirectLightPowerPerArea = 0.5f*vec3( 1, 1, 1 );

// Normal and reflectance values to be passed to PixelShader
vec4 currentNormal;
vec3 currentReflectance;

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
  LoadTestModel( triangles );
  //std::vector<glm::ivec2> points;

  /* Clear buffer */
  memset(screen->buffer, 0, screen->height*screen->width*sizeof(uint32_t));

  // Clear depth buffer.
  for (int y=0; y<SCREEN_HEIGHT; ++y) {
    for (int x=0; x<SCREEN_WIDTH; ++x) {
      depthBuffer[y][x] = 0;
    }
  }

  // Go through all triangles
  for( uint32_t i=0; i<triangles.size(); ++i ) {

     //points.clear();
     vector<Vertex> vertices(3);

     vertices[0].position = triangles[i].v0;
     vertices[1].position = triangles[i].v1;
     vertices[2].position = triangles[i].v2;

     // Set current colour, reflectance and normal to the colour of the triangle.
     currentColor = triangles[i].color;
     currentNormal = triangles[i].normal;
     //currentReflectance = ??;

     // Draw and colour the triangles.
     DrawPolygon( screen, vertices );
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
          lightPos += vec4(0, 0.1, 0, 0);
        break;
        case SDLK_e:
          lightPos += vec4(0, -0.1, 0, 0);
        break;
        // CAMERA MOVEMENT
        case SDLK_UP:
        cameraPos += vec4(0, 0, 0.1, 0);
        /* Move camera forward */
        break;
        case SDLK_DOWN:
        cameraPos += vec4(0, 0, -0.1, 0);
        /* Move camera backwards */
        break;
        case SDLK_LEFT:
        cameraPos += vec4(-0.1, 0, 0, 0);
        /* Move camera left */
        break;
        case SDLK_RIGHT:
        cameraPos += vec4(0.1, 0, 0, 0);
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
    if (i == 2) b = vertexPixels[0];
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
  float x3D = v.position.x - cameraPos[0];
  float y3D = v.position.y - cameraPos[1];
  float z3D = v.position.z - cameraPos[2];


  float x = (focalLength * (x3D / z3D)) + (SCREEN_WIDTH/2);
  float y = (focalLength * (y3D / z3D)) + (SCREEN_HEIGHT/2);

  int x2D = static_cast<int>(x);
  int y2D = static_cast<int>(y);

  p.x = x2D;
  p.y = y2D;
  p.zinv = 1/z3D;

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

     result[i].pos3d.x = a.pos3d.x + (stepPos3dX * i);
     result[i].pos3d.y = a.pos3d.y + (stepPos3dY * i);
     result[i].pos3d.z = a.pos3d.z + (stepPos3dZ * i);
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
  if (p.zinv > depthBuffer[y][x]) {
    depthBuffer[y][x] = p.zinv;
    PutPixelSDL( screen, x, y, currentColor * calculateIllumination(p, currentNormal) );
    // printf("%f, %f, %f\n", p.illumination.x, p.illumination.y, p.illumination.z );
  }
}

vec3 calculateIllumination(const Pixel& p, vec4 currentNormal) {
  vec4 r = lightPos - p.pos3d;
  vec3 vec3r = vec3(r);
  float r_magnitude = pow(r[0],2) + pow(r[1],2) + pow(r[2],2);

  vec3 normal = vec3(currentNormal);

  float vecProduct = glm::dot(vec3r, normal);
  vec3 D = (lightPower * glm::max(vecProduct, 0.0f)) / (float)(4.0f * M_PI * r_magnitude);

  vec3 totalIllumination = D + indirectLightPowerPerArea;

  return totalIllumination;

}
