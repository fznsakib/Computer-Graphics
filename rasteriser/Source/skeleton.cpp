#include <iostream>
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

/* ----------------------------------------------------------------------------*/
/* FUNCTIONS                                                                   */
/* ----------------------------------------------------------------------------*/

bool Update();
void Draw(screen* screen);
void VertexShader( const vec4& v, glm::ivec2& p );
void Interpolate( glm::ivec2 a, glm::ivec2 b, vector<glm::ivec2>& result );
void DrawLineSDL( screen* screen, vector<glm::ivec2> line, vec3 color );
void DrawPolygon( screen* screen, const vector<vec4>& vertices);
void ComputePolygonRows( const vector<glm::ivec2>& vertexPixels, vector<glm::ivec2>& leftPixels, vector<glm::ivec2>& rightPixels);
void DrawPolygonRows( screen* screen, const vector<glm::ivec2>& leftPixels, const vector<glm::ivec2>& rightPixels );


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

  /* TO TEST INTERPOLATION VALUES */
  // vector<glm::ivec2> vertexPixels(3);
  // vertexPixels[0] = glm::ivec2(10, 5);
  // vertexPixels[1] = glm::ivec2( 5,10);
  // vertexPixels[2] = glm::ivec2(15,15);
  // vector<glm::ivec2> leftPixels;
  // vector<glm::ivec2> rightPixels;
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
  std::vector<glm::ivec2> points;

  /* Clear buffer */
  memset(screen->buffer, 0, screen->height*screen->width*sizeof(uint32_t));

  // Go through all triangles
  for( uint32_t i=0; i<triangles.size(); ++i ) {

     points.clear();
     vector<vec4> vertices(3);

     vertices[0] = triangles[i].v0;
     vertices[1] = triangles[i].v1;
     vertices[2] = triangles[i].v2;

     // Set current colour to the colour of the triangle.
     currentColor = triangles[i].color;

     // Draw and colour the triangles.
     DrawPolygon( screen, vertices );

    //  // Draw vertices
    // for(int v=0; v<3; ++v) {
    //    glm::ivec2 projPos;
    //    VertexShader( vertices[v], projPos );
    //    vec3 color(1,1,1);
    //    PutPixelSDL( screen, projPos.x, projPos.y, color );
    //
    //    points.push_back(projPos);
    // }
    //
    // glm::ivec2 a;
    // glm::ivec2 b;
    // // Get edges between vertices
    // for(int i = 0; i < 3; i++) {
    //   a = points[i];
    //   if (i == 2) b = points[0];
    //   else b = points[i + 1];
    //
    //   // Calculate number of pixels to draw in interpolation
    //   glm::ivec2 delta = glm::abs( a - b );
    //   int pixels = glm::max( delta.x, delta.y ) + 1;
    //
    //   // Populate vector with pixels to be drawn
    //   vector<glm::ivec2> line( pixels );
    //   Interpolate( a, b, line );
    //
    //   // Draw all pixels within line
    //   DrawLineSDL(screen, line, vec3(1,1,1));
    // }
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
        /* --------------- ADDED KEYPRESS CAMERA MOVEMENTS -------------------*/
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
void DrawPolygon(screen* screen, const vector<vec4>& vertices) {
  int V = vertices.size();

  vector<glm::ivec2> vertexPixels(V);
  for (int i=0; i<V; ++i) {
    VertexShader( vertices[i], vertexPixels[i] );
  }
  vector<glm::ivec2> leftPixels;
  vector<glm::ivec2> rightPixels;
  ComputePolygonRows( vertexPixels, leftPixels, rightPixels );
  DrawPolygonRows( screen, leftPixels, rightPixels );
}

void ComputePolygonRows( const vector<glm::ivec2>& vertexPixels, vector<glm::ivec2>& leftPixels, vector<glm::ivec2>& rightPixels) {
  int max = -numeric_limits<int>::max();
  int min = +numeric_limits<int>::max();
  int rows = 0;

  // Find min and max y-value of polygon and compute number of rows it occupies.
  for (int i = 0; i < vertexPixels.size(); ++i) {
    glm::ivec2 currentVertex = vertexPixels[i];
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
  glm::ivec2 a;
  glm::ivec2 b;
  for(int i = 0; i < vertexPixels.size(); i++) {
    a = vertexPixels[i];
    if (i == 2) b = vertexPixels[0];
    else b = vertexPixels[i + 1];

    // Calculate number of pixels to draw in interpolation
    glm::ivec2 delta = glm::abs( a - b );
    int pixels = glm::max( delta.x, delta.y ) + 1;

    // Populate vector with pixels to be drawn
    vector<glm::ivec2> line( pixels );
    Interpolate( a, b, line );

    /* interpolation changed to floor() instead of round() to achieve correct vals. */
    /*----------------------------------------------*/
    for (int j = 0; j < line.size(); ++j) {
      if (line[j].x <= leftPixels[line[j].y - min].x && line[j].y - min > 0) { // SEG FAULT FIXED HERE
        leftPixels[line[j].y - min].x = line[j].x;
      }
      if (line[j].x >= rightPixels[line[j].y - min].x && line[j].y - min > 0) { // SEG FAULT FIXED HERE
        rightPixels[line[j].y - min].x = line[j].x;
      }
    }
    /*----------------------------------------------*/
  }
}

void DrawPolygonRows( screen* screen, const vector<glm::ivec2>& leftPixels, const vector<glm::ivec2>& rightPixels ) {
  for (int i = 0; i < leftPixels.size(); ++i) {
    for (int j = leftPixels[i].x; j < rightPixels[i].x; ++j) {
      PutPixelSDL(screen, j, leftPixels[i].y, currentColor);
    }
  }
}

void VertexShader( const vec4& v, glm::ivec2& p ) {
  float x3D = v[0] - cameraPos[0];
  float y3D = v[1] - cameraPos[1];
  float z3D = v[2] - cameraPos[2];


  float x = (focalLength * (x3D / z3D)) + (SCREEN_WIDTH/2);
  float y = (focalLength * (y3D / z3D)) + (SCREEN_HEIGHT/2);

  int x2D = static_cast<int>(x);
  int y2D = static_cast<int>(y);

  p[0] = x2D;
  p[1] = y2D;

  // std::cout << x2D << " and " << y2D << '\n';
}

void Interpolate( glm::ivec2 a, glm::ivec2 b, vector<glm::ivec2>& result ) {
  int N = result.size();
  glm::vec2 step = glm::vec2(b-a) / float(max(N-1,1));
  glm::vec2 current( a );

  for( int i=0; i<N; ++i ) {
     result[i] = floor(current); // Change to floor to achieve same interpolation result as notes.
     current += step;
  }
}

void DrawLineSDL( screen* screen, vector<glm::ivec2> line, vec3 color ) {
    for(int i = 0; i < line.size(); i++) {
      PutPixelSDL(screen, line[i].x, line[i].y, color);
    }
}
