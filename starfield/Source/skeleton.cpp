#include <iostream>
#include <glm/glm.hpp>
#include <SDL.h>
#include "SDLauxiliary.h"
#include "TestModel.h"
#include <stdint.h>

using namespace std;
using glm::vec3;
using glm::mat3;

#define SCREEN_WIDTH 320
#define SCREEN_HEIGHT 256
#define FULLSCREEN_MODE true


/* ----------------------------------------------------------------------------*/
/* GLOBAL VARIABLES                                                            */
/* ----------------------------------------------------------------------------*/

int t;

/* ----------------------------------------------------------------------------*/
/* FUNCTIONS                                                                   */
/* ----------------------------------------------------------------------------*/

void Update();
void Draw(screen* screen);
void Interpolate( float a, float b, vector<float>& result );
void Interpolate( vec3 a, vec3 b, vector<vec3>& result );

int main( int argc, char* argv[] ) {

  screen *screen = InitializeSDL( SCREEN_WIDTH, SCREEN_HEIGHT, FULLSCREEN_MODE );
  t = SDL_GetTicks();	/*Set start value for timer.*/

  while( NoQuitMessageSDL() ) {
    Draw(screen);
    Update();
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

  vec3 topLeft(1,0,0); //red
  vec3 topRight(0,0,1); // blue
  vec3 bottomRight(0,1,0); // green
  vec3 bottomLeft(1,1,0); // yellow

  vector<vec3> leftSide( SCREEN_HEIGHT );
  vector<vec3> rightSide( SCREEN_HEIGHT );
  Interpolate( topLeft, bottomLeft, leftSide );
  Interpolate( topRight, bottomRight, rightSide );

  // vec3 colour(1.0,0.0,0.0);
  for(int y = 0; y < SCREEN_HEIGHT; y++) {
    vector<vec3> row( SCREEN_WIDTH );
    Interpolate(leftSide[y], rightSide[y], row);

    for(int x = 0; x < SCREEN_WIDTH; x++) {
      PutPixelSDL(screen, x, y, row[x]);
    }
  }
}

/*Place updates of parameters here*/
void Update() {
  /* Compute frame time */
  int t2 = SDL_GetTicks();
  //float dt = float(t2-t);
  t = t2;
  /*Good idea to remove this*/
  // std::cout << "Render time: " << dt << " ms." << std::endl;
  /* Update variables*/
}

void Interpolate( float a, float b, vector<float>& result ) {

  if (result.size() == 1) {
    result[0] = (a + b)/2;
  }
  else {
    float section = (b - a)/(result.size() - 1);

    for (int i = 0; i < result.size(); i++) {
      result[i] = a + (section * i);
    }
  }
}

void Interpolate( vec3 a, vec3 b, vector<vec3>& result ) {
  if (result.size() == 1) {
    result[0].x = (a.x + b.x)/2;
    result[0].y = (a.y + b.y)/2;
    result[0].z = (a.z + b.z)/2;
  }

  float sectionX = (b.x - a.x)/(result.size() - 1);
  float sectionY = (b.y - a.y)/(result.size() - 1);
  float sectionZ = (b.z - a.z)/(result.size() - 1);

  for (int i = 0; i < result.size(); i++) {
    result[i].x = a.x + (sectionX * i);
    result[i].y = a.y + (sectionY * i);
    result[i].z = a.z + (sectionZ * i);
  }
}
