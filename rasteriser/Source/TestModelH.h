#ifndef TEST_MODEL_CORNEL_BOX_H
#define TEST_MODEL_CORNEL_BOX_H

// Defines a simple test model: The Cornel Box

#include <glm/glm.hpp>
#include <vector>

int setting = 3;

// Used to describe a triangular surface:
class Triangle
{
public:
	glm::vec4 v0;
	glm::vec4 v1;
	glm::vec4 v2;
	glm::vec4 normal;
	glm::vec3 color;
	// Set textures with this value. 1 = marble, 2 = metal grill, 3 = woven
	int texture = 0;
	// Index to object. 0 = back, 1 = ceiling, 2 = floor, 3 = leftwall, 4 = rightwall
	int index;

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
	Triangle floor1( C, B, A, green );
	floor1.texture = setting;
	floor1.index = 2;
	room.push_back( floor1 );
	Triangle floor2( C, D, B, green );
	floor2.texture = setting;
	floor2.index = 2;
	room.push_back( floor2 );

	// Left wall
	Triangle leftwall1 ( A, E, C, purple );
	leftwall1.texture = setting;
	leftwall1.index = 3;
	room.push_back( leftwall1 );
	Triangle leftwall2 ( C, E, G, purple );
	leftwall2.texture = setting;
	leftwall2.index = 3;
	room.push_back( leftwall2 );

	// Right wall
	Triangle rightwall1( F, B, D, yellow );
	rightwall1.texture = setting;
	rightwall1.index = 4;
	room.push_back( rightwall1  );
	Triangle rightwall2( H, F, D, yellow );
	rightwall2.texture = setting;
	rightwall2.index = 4;
	room.push_back( rightwall2 );

	// Ceiling
	Triangle ceiling1 ( E, F, G, cyan );
	ceiling1.texture = setting;
	ceiling1.index = 1;
	room.push_back( ceiling1 );
	Triangle ceiling2( F, H, G, cyan );
	ceiling2.texture = setting;
	ceiling2.index = 1;
	room.push_back( ceiling2 );

	// Back wall
	Triangle backwall1( G, D, C, vec3(0.03529f, 0.7843f, 0.8078f) );
	backwall1.texture = setting;
	backwall1.index = 0;
	room.push_back( backwall1 );
	Triangle backwall2( G, H, D, vec3(0.03529f, 0.7843f, 0.8078f) );
	backwall2.texture = setting;
	backwall2.index = 0;
	room.push_back( backwall2 );

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
	Triangle a(E,B,A,red);
	a.texture = 1;
	a.index = 0;
	boxes.push_back( a );
	Triangle s(E,F,B,red);
	s.texture = 1;
	s.index = 0;
	boxes.push_back( s );

	// Front
	Triangle d(F,D,B,red);
	d.texture = 1;
	d.index = 4;
	boxes.push_back( d );
	Triangle f(F,H,D,red);
	f.texture = 1;
	f.index = 4;
	boxes.push_back( f );

	// BACK
	Triangle g(H,C,D,red);
	g.texture = 1;
	g.index = 0;
	boxes.push_back( g );
	Triangle h(H,G,C,red);
	h.texture = 1;
	h.index = 0;
	boxes.push_back( h );

	// LEFT
	Triangle j(G,E,C,red);
	j.texture = 1;
	j.index = 3;
	boxes.push_back( j );
	Triangle k(E,A,C,red);
	k.texture = 1;
	k.index = 3;
	boxes.push_back( k );

	// TOP
	Triangle l(G,F,E,red);
	l.texture = 1;
	l.index = 1;
	boxes.push_back( l );
	Triangle q(G,H,F,red);
	q.texture = 1;
	q.index = 1;
	boxes.push_back( q );

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
	Triangle front_tallBlock1(E,B,A,blue);
	front_tallBlock1.texture = 3;
	front_tallBlock1.index = 0;
	boxes.push_back( front_tallBlock1 );
	Triangle front_tallBlock2(E,F,B,blue);
	front_tallBlock2.texture = 3;
	front_tallBlock2.index = 0;
	boxes.push_back( front_tallBlock2 );

	// Front
	Triangle bottom_tallBlock1(F,D,B,blue);
	bottom_tallBlock1.texture = 3;
	bottom_tallBlock1.index = 4;
	boxes.push_back( bottom_tallBlock1 );
	Triangle bottom_tallBlock2(F,H,D,blue);
	bottom_tallBlock2.texture = 3;
	bottom_tallBlock2.index = 4;
	boxes.push_back( bottom_tallBlock2 );

	// BACK
	Triangle back_tallBlock1(H,C,D,blue);
	back_tallBlock1.texture = 3;
	back_tallBlock1.index = 0;
	boxes.push_back( back_tallBlock1 );
	Triangle back_tallBlock2(H,G,C,blue);
	back_tallBlock2.texture = 3;
	back_tallBlock2.index = 0;
	boxes.push_back( back_tallBlock2 );

	// LEFT
	Triangle left_tallBlock1(G,E,C,blue);
	left_tallBlock1.texture = 3;
	left_tallBlock1.index = 3;
	boxes.push_back( left_tallBlock1 );
	Triangle left_tallBlock2(E,A,C,blue);
	left_tallBlock2.texture = 3;
	left_tallBlock2.index = 3;
	boxes.push_back( left_tallBlock2 );

	// TOP
	Triangle top_tallBlock1(G,F,E,blue);
	top_tallBlock1.texture = 3;
	top_tallBlock1.index = 1;
	boxes.push_back( top_tallBlock1 );
	Triangle top_tallBlock2(G,H,F,blue);
	top_tallBlock2.texture = 3;
	top_tallBlock1.index = 1;
	boxes.push_back( top_tallBlock2 );


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
