#ifndef TEST_MODEL_CORNEL_BOX_H
#define TEST_MODEL_CORNEL_BOX_H

// Defines a simple test model: The Cornel Box

#include <glm/glm.hpp>
#include <vector>

// Used to describe a triangular surface:
class Triangle
{
public:
	glm::vec4 v0;
	glm::vec4 v1;
	glm::vec4 v2;
	glm::vec4 normal;
	glm::vec3 color;
	bool mirror = false;

	Triangle( glm::vec4 v0, glm::vec4 v1, glm::vec4 v2, glm::vec3 color )
		: v0(v0), v1(v1), v2(v2), color(color)
	{
		ComputeNormal();
	}

	void ComputeNormal()
	{
	  glm::vec3 e1 = glm::vec3(v1.x-v0.x,v1.y-v0.y,v1.z-v0.z);
	  glm::vec3 e2 = glm::vec3(v2.x-v0.x,v2.y-v0.y,v2.z-v0.z);
	  glm::vec3 normal3 = glm::normalize( glm::cross( e2, e1 ) );
	  normal.x = normal3.x;
	  normal.y = normal3.y;
	  normal.z = normal3.z;
	  normal.w = 1.0;
	}
};

// Loads the Cornell Box. It is scaled to fill the volume:
// -1 <= x <= +1
// -1 <= y <= +1
// -1 <= z <= +1
void LoadTestModel( std::vector<Triangle>& room, std::vector<Triangle>& boxes )
{
	using glm::vec3;
	using glm::vec4;

	// Defines colors:
	vec3 red(    0.75f, 0.15f, 0.15f );
	vec3 yellow( 0.75f, 0.75f, 0.15f );
	vec3 green(  0.15f, 0.75f, 0.15f );
	vec3 cyan(   0.15f, 0.75f, 0.75f );
	vec3 blue(   0.15f, 0.15f, 0.75f );
	vec3 purple( 0.75f, 0.15f, 0.75f );
	vec3 white(  0.75f, 0.75f, 0.75f );
	vec3 wood( 10.0f, 10.0f, 10.0f);

	room.clear();
	boxes.clear();
	room.reserve( 5*2 );
	boxes.reserve( 20 );


	// ---------------------------------------------------------------------------
	// Room

	float L = 555;			// Length of Cornell Box side.

	vec4 A(L,0,0,1);
	vec4 B(0,0,0,1);
	vec4 C(L,0,L,1);
	vec4 D(0,0,L,1);

	vec4 E(L,L,0,1);
	vec4 F(0,L,0,1);
	vec4 G(L,L,L,1);
	vec4 H(0,L,L,1);

	// Floor:
	room.push_back( Triangle( C, B, A, green ) );
	room.push_back( Triangle( C, D, B, green ) );

	// Left wall
	Triangle leftwall1 ( A, E, C, purple );
	leftwall1.mirror = true;
	room.push_back( leftwall1 );
	Triangle leftwall2 ( C, E, G, purple );
	leftwall2.mirror = true;
	room.push_back( leftwall2 );

	// Right wall
	room.push_back( Triangle( F, B, D, yellow ) );
	room.push_back( Triangle( H, F, D, yellow ) );

	// Ceiling
	Triangle ceiling1 ( E, F, G, cyan );
	ceiling1.mirror = true;
	room.push_back( ceiling1 );
	Triangle ceiling2( F, H, G, cyan );
	ceiling2.mirror = true;
	room.push_back( ceiling2 );

	// Back wall
	room.push_back( Triangle( G, D, C, white ) );
	room.push_back( Triangle( G, H, D, white ) );

	// ---------------------------------------------------------------------------
	// Short block

	A = vec4(290,0,114,1);
	B = vec4(130,0, 65,1);
	C = vec4(240,0,272,1);
	D = vec4( 82,0,225,1);

	E = vec4(290,165,114,1);
	F = vec4(130,165, 65,1);
	G = vec4(240,165,272,1);
	H = vec4( 82,165,225,1);

	// Front
	boxes.push_back( Triangle(E,B,A,red) );
	boxes.push_back( Triangle(E,F,B,red) );

	// Front
	boxes.push_back( Triangle(F,D,B,red) );
	boxes.push_back( Triangle(F,H,D,red) );

	// BACK
	boxes.push_back( Triangle(H,C,D,red) );
	boxes.push_back( Triangle(H,G,C,red) );

	// LEFT
	boxes.push_back( Triangle(G,E,C,red) );
	boxes.push_back( Triangle(E,A,C,red) );

	// TOP
	boxes.push_back( Triangle(G,F,E,red) );
	boxes.push_back( Triangle(G,H,F,red) );

	// ---------------------------------------------------------------------------
	// Tall block

	A = vec4(423,0,247,1);
	B = vec4(265,0,296,1);
	C = vec4(472,0,406,1);
	D = vec4(314,0,456,1);

	E = vec4(423,330,247,1);
	F = vec4(265,330,296,1);
	G = vec4(472,330,406,1);
	H = vec4(314,330,456,1);

	// Front
	boxes.push_back( Triangle(E,B,A,blue) );
	boxes.push_back( Triangle(E,F,B,blue) );

	// Front
	boxes.push_back( Triangle(F,D,B,blue) );
	boxes.push_back( Triangle(F,H,D,blue) );

	// BACK
	boxes.push_back( Triangle(H,C,D,blue) );
	boxes.push_back( Triangle(H,G,C,blue) );

	// LEFT
	boxes.push_back( Triangle(G,E,C,blue) );
	boxes.push_back( Triangle(E,A,C,blue) );

	// TOP
	boxes.push_back( Triangle(G,F,E,blue) );
	boxes.push_back( Triangle(G,H,F,blue) );


	// ----------------------------------------------
	// Scale to the volume [-1,1]^3

	for( size_t i=0; i<room.size(); ++i )
	{
		room[i].v0 *= 2/L;
		room[i].v1 *= 2/L;
		room[i].v2 *= 2/L;

		room[i].v0 -= vec4(1,1,1,1);
		room[i].v1 -= vec4(1,1,1,1);
		room[i].v2 -= vec4(1,1,1,1);

		room[i].v0.x *= -1;
		room[i].v1.x *= -1;
		room[i].v2.x *= -1;

		room[i].v0.y *= -1;
		room[i].v1.y *= -1;
		room[i].v2.y *= -1;

		room[i].v0.w = 1.0;
		room[i].v1.w = 1.0;
		room[i].v2.w = 1.0;

		room[i].ComputeNormal();
	}

	for( size_t i=0; i<boxes.size(); ++i )
	{
		boxes[i].v0 *= 2/L;
		boxes[i].v1 *= 2/L;
		boxes[i].v2 *= 2/L;

		boxes[i].v0 -= vec4(1,1,1,1);
		boxes[i].v1 -= vec4(1,1,1,1);
		boxes[i].v2 -= vec4(1,1,1,1);

		boxes[i].v0.x *= -1;
		boxes[i].v1.x *= -1;
		boxes[i].v2.x *= -1;

		boxes[i].v0.y *= -1;
		boxes[i].v1.y *= -1;
		boxes[i].v2.y *= -1;

		boxes[i].v0.w = 1.0;
		boxes[i].v1.w = 1.0;
		boxes[i].v2.w = 1.0;

		boxes[i].ComputeNormal();
	}
}

#endif
