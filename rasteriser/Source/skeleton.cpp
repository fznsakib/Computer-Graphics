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


/* ----------------------------------------------------------------------------*/
/* FUNCTIONS                                                                   */
/* ----------------------------------------------------------------------------*/

bool Update();
void Draw(screen* screen);
void VertexShader( const vec4& v, glm::ivec2& p );

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
void Draw(screen* screen)
{
  vector<Triangle> triangles;
  LoadTestModel( triangles );

  /* Clear buffer */
  memset(screen->buffer, 0, screen->height*screen->width*sizeof(uint32_t));
  for( uint32_t i=0; i<triangles.size(); ++i ) {
     vector<vec4> vertices(3);

     vertices[0] = triangles[i].v0;
     vertices[1] = triangles[i].v1;
     vertices[2] = triangles[i].v2;

     for(int v=0; v<3; ++v) {
       glm::ivec2 projPos;
       VertexShader( vertices[v], projPos );
       vec3 color(1,1,1);
       PutPixelSDL( screen, projPos.x, projPos.y, color );
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

  std::cout << x2D << " and " << y2D << '\n';
}
