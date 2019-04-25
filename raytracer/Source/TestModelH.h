#ifndef TEST_MODEL_CORNEL_BOX_H
#define TEST_MODEL_CORNEL_BOX_H

// Defines a simple test model: The Cornel Box

#include <glm/glm.hpp>
#include <vector>

using glm::vec3;
using glm::vec4;


// Used to describe a spherical surface:
class Sphere
{
public:
		Sphere(const vec3 &c, const float &r, const vec3 &color) :
           radius(r), radiusSquared(r * r), centre(c), color(color) {}

		float radius, radiusSquared;
		vec3 centre, color, normal;
    // vec4 normal;

		bool solveQuadratic(const float &a, const float &b, const float &c, float &x0, float &x1) const
		{
				// Check discriminant first to see how many solutions there are
		    float discriminant = (b * b) - (4 * a * c);
		    if (discriminant < 0) return false;
		    else if (discriminant == 0) x0 = x1 = - 0.5 * b / a;
		    else {
					float q;
					if (b > 0) q = -0.5 * (b + sqrt(discriminant));
					else 			 q = -0.5 * (b - sqrt(discriminant));
	        x0 = q / a;
	        x1 = c / q;
		    }
		    if (x0 > x1) std::swap(x0, x1);

		    return true;
		}

		// To compute the intersection of ray with a sphere
    bool intersect(const vec3 &start, const vec3 &dir, float &t) const
    {
        float t0, t1; // solutions for t if the ray intersects

				// Find variables to get solutions from quadratic equation
				vec3 L = start - centre;
        float a = glm::dot(dir, dir);
        float b = 2 * glm::dot(dir, L);
        float c = glm::dot(L, L) - radiusSquared;
        if (!solveQuadratic(a, b, c, t0, t1)) return false;

				if (t0 > t1) std::swap(t0, t1);

				// If t0 is negative, there is only one root. Use t0 only then
        if (t0 < 0) {
            t0 = t1;
						// If both t0 and t1 are negative, then there is no intersection
            if (t0 < 0) return false;
        }

        t = t0;

        return true;
    }

    void getNormal(const vec3 &position, vec3 &normal) const
    {
        // In this particular case, the normal is similar to a point on a unit sphere
        // centred around the origin. We can thus use the normal coordinates to compute
        // the spherical coordinates of where the ray hits the surface of the sphere
        normal = position - centre;
        normal = glm::normalize(normal);
    }

};

// Used to describe a triangular surface:
class Triangle
{
public:

	Triangle( glm::vec4 v0, glm::vec4 v1, glm::vec4 v2, glm::vec3 color )
		: v0(v0), v1(v1), v2(v2), color(color)
	{
		ComputeNormal();
	}

	glm::vec4 v0;
	glm::vec4 v1;
	glm::vec4 v2;
	glm::vec4 normal;
	glm::vec3 color;

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


	// To compute the intersection of ray with a triangle
	bool intersect(const vec3 &start, const vec3 &dir, float &t) const
  {
      // code to compute the intersection of a ray with triangle
			if (t< 0) return false;
      else return true;
  }
};

// Loads the Cornell Box. It is scaled to fill the volume:
// -1 <= x <= +1
// -1 <= y <= +1
// -1 <= z <= +1
void LoadTestModel( std::vector<Triangle>& triangles, std::vector<Sphere>& spheres )
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

	triangles.clear();
	triangles.reserve( 5*2*3 );

	// ---------------------------------------------------------------------------
	// Room
	// ---------------------------------------------------------------------------

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
	triangles.push_back( Triangle( C, B, A, green ) );
	triangles.push_back( Triangle( C, D, B, green ) );

	// Left wall
	triangles.push_back( Triangle( A, E, C, purple ) );
	triangles.push_back( Triangle( C, E, G, purple ) );

	// Right wall
	triangles.push_back( Triangle( F, B, D, yellow ) );
	triangles.push_back( Triangle( H, F, D, yellow ) );

	// Ceiling
	triangles.push_back( Triangle( E, F, G, cyan ) );
	triangles.push_back( Triangle( F, H, G, cyan ) );

	// Back wall
	triangles.push_back( Triangle( G, D, C, white ) );
	triangles.push_back( Triangle( G, H, D, white ) );

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
	triangles.push_back( Triangle(E,B,A,red) );
	triangles.push_back( Triangle(E,F,B,red) );

	// Front
	triangles.push_back( Triangle(F,D,B,red) );
	triangles.push_back( Triangle(F,H,D,red) );

	// BACK
	triangles.push_back( Triangle(H,C,D,red) );
	triangles.push_back( Triangle(H,G,C,red) );

	// LEFT
	triangles.push_back( Triangle(G,E,C,red) );
	triangles.push_back( Triangle(E,A,C,red) );

	// TOP
	triangles.push_back( Triangle(G,F,E,red) );
	triangles.push_back( Triangle(G,H,F,red) );


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
	triangles.push_back( Triangle(E,B,A,blue) );
	triangles.push_back( Triangle(E,F,B,blue) );

	// Front
	triangles.push_back( Triangle(F,D,B,blue) );
	triangles.push_back( Triangle(F,H,D,blue) );

	// BACK
	// triangles.push_back( Triangle(H,C,D,blue) );
	// triangles.push_back( Triangle(H,G,C,blue) );

	// LEFT
	triangles.push_back( Triangle(G,E,C,blue) );
	triangles.push_back( Triangle(E,A,C,blue) );

	// TOP
	triangles.push_back( Triangle(G,F,E,blue) );
	triangles.push_back( Triangle(G,H,F,blue) );


	// ----------------------------------------------
	// Scale to the volume [-1,1]^3

	for( size_t i=0; i<triangles.size(); ++i )
	{
		triangles[i].v0 *= 2/L;
		triangles[i].v1 *= 2/L;
		triangles[i].v2 *= 2/L;

		triangles[i].v0 -= vec4(1,1,1,1);
		triangles[i].v1 -= vec4(1,1,1,1);
		triangles[i].v2 -= vec4(1,1,1,1);

		triangles[i].v0.x *= -1;
		triangles[i].v1.x *= -1;
		triangles[i].v2.x *= -1;

		triangles[i].v0.y *= -1;
		triangles[i].v1.y *= -1;
		triangles[i].v2.y *= -1;

		triangles[i].v0.w = 1.0;
		triangles[i].v1.w = 1.0;
		triangles[i].v2.w = 1.0;

		triangles[i].ComputeNormal();
	}


	// ----------------------------------------------
	// Spheres

	vec3 centre(-0.45, 0.6, -0.6);
  float radius = 0.3;
	spheres.push_back(Sphere(centre, radius, red));

}

#endif
