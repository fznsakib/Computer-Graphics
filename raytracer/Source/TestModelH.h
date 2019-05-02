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
		Sphere(const vec3 &c, const float &r, const vec3 &color, const vec3 &material) :
           radius(r), radiusSquared(r * r), centre(c), color(color), material(material) {}

		float radius, radiusSquared;
		vec3 centre, color, normal, material;
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

    void getNormal(const vec3 &position, vec4 &normal) const
    {
        // In this particular case, the normal is similar to a point on a unit sphere
        // centred around the origin. We can thus use the normal coordinates to compute
        // the spherical coordinates of where the ray hits the surface of the sphere
				vec3 normal3 = position - centre;
				normal3 = glm::normalize(normal3);

        normal = vec4(normal3.x, normal3.y, normal3.z, 1.0f);
    }

};

// Used to describe a triangular surface:
class Triangle
{
public:

	Triangle( glm::vec4 v0, glm::vec4 v1, glm::vec4 v2, glm::vec3 color,  glm::vec3 material )
		: v0(v0), v1(v1), v2(v2), color(color), material(material)
	{
		ComputeNormal();
	}

	glm::vec4 v0;
	glm::vec4 v1;
	glm::vec4 v2;
	glm::vec4 normal;
	glm::vec3 color;
	glm::vec3 material;

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
	bool intersect(const vec3 &start, const vec3 &dir, float &t, const float &closestIntersectionDistance) const
  {
		vec4 start4 = vec4(start[0], start[1], start[2], 0);
    vec3 dir3 = vec3(dir[0], dir[1], dir[2]);
		float bound = std::numeric_limits<float>::max();


		vec3 e1 = vec3(v1.x - v0.x ,v1.y - v0.y,v1.z - v0.z);
		vec3 e2 = vec3(v2.x - v0.x ,v2.y - v0.y,v2.z - v0.z);
		// vec3 b = vec3(start.x-v0.x, start.y-v0.y, start.z-v0.z);

		// Currently: 55000ms without -O3, 600ms with -O3
		// With Cramer's rule: 47000ms without -O3, 430ms with -O3
		glm::mat3 system( -dir3, e1, e2 );

		// Use Cramer's rule to produce just t (x[0]). If t is negative, then
		// return false. Else, continue with calculation

		// Using the linear equation:
		// (-dir, e1, e2) (t, u, v) = s - v0
		vec4 solutions = start4 - v0;
		vec3 solutions3 = vec3(solutions[0], solutions[1], solutions[2]);

		// Solving t in using Cramer's rule:
		// -dir[0]t + e1[0]u + e2[0]v = solutions[0]
		// -dir[1]t + e1[1]u + e2[1]v = solutions[1]
		// -dir[2]t + e1[2]u + e2[2]v = solutions[2]

		// Where t = det(system_t)/det(system)
		glm::mat3 system_t( solutions3, e1, e2 );
		t = glm::determinant(system_t)/glm::determinant(system);
		// float distance = t * glm::length(dir3);


		// If the distance is negative, the intersection does not occur, so carry on to the next one.
		if (t < 0.0f) return false;
		// Do distance checks earlier to move on quicker if possible
		else if (t >= closestIntersectionDistance || t > bound) return false;


		// Similar to t, calculate u and v
		glm::mat3 system_u( -dir3, solutions3, e2 );
		float u = glm::determinant(system_u)/glm::determinant(system);

		glm::mat3 system_v( -dir3, e1, solutions3 );
		float v = glm::determinant(system_v)/glm::determinant(system);

		// x = (t, u, v)
		// vec3 x = glm::inverse( system ) * b;

		bool check = (u >= 0) && (v >= 0) && ((u + v) <= 1);

		if (check) return true;
		else return false;
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
	vec3 black(  0.0f, 0.0f, 0.0f );

	vec3 beige(0.85f, 0.85f, 0.7f);

	// Define material:
	// Material = (diffuse, specular, transmission)
	vec3 matte(0.5f, 0.0f, 0.0f);
	vec3 chrome(0.0f, 0.8f, 0.0f);
	vec3 glass(0.0f, 0.1f, 0.8f);


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
	triangles.push_back( Triangle( C, B, A, white, matte ) );
	triangles.push_back( Triangle( C, D, B, white, matte ) );

	// Left wall
	triangles.push_back( Triangle( A, E, C, cyan, matte ) );
	triangles.push_back( Triangle( C, E, G, cyan, matte ) );

	// Right wall
	triangles.push_back( Triangle( F, B, D, cyan, matte ) );
	triangles.push_back( Triangle( H, F, D, cyan, matte ) );

	// Ceiling
	triangles.push_back( Triangle( E, F, G, white, matte ) );
	triangles.push_back( Triangle( F, H, G, white, matte ) );

	// Back wall
	triangles.push_back( Triangle( G, D, C, red, matte ) );
	triangles.push_back( Triangle( G, H, D, red, matte ) );

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

	// // Front
	// triangles.push_back( Triangle(E,B,A,red, matte) );
	// triangles.push_back( Triangle(E,F,B,red, matte) );
	//
	// // Front
	// triangles.push_back( Triangle(F,D,B,red, matte) );
	// triangles.push_back( Triangle(F,H,D,red, matte) );
	//
	// // BACK
	// triangles.push_back( Triangle(H,C,D,red, matte) );
	// triangles.push_back( Triangle(H,G,C,red, matte) );
	//
	// // LEFT
	// triangles.push_back( Triangle(G,E,C,red, matte) );
	// triangles.push_back( Triangle(E,A,C,red, matte) );
	//
	// // TOP
	// triangles.push_back( Triangle(G,F,E,red, matte) );
	// triangles.push_back( Triangle(G,H,F,red, matte) );


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
	// triangles.push_back( Triangle(E,B,A,blue, matte) );
	// triangles.push_back( Triangle(E,F,B,blue, matte) );

	// Front
	// triangles.push_back( Triangle(F,D,B,blue, matte) );
	// triangles.push_back( Triangle(F,H,D,blue, matte) );

	// BACK
	// triangles.push_back( Triangle(H,C,D,blue) );
	// triangles.push_back( Triangle(H,G,C,blue) );

	// LEFT
	// triangles.push_back( Triangle(G,E,C,blue, matte) );
	// triangles.push_back( Triangle(E,A,C,blue, matte) );

	// TOP
	// triangles.push_back( Triangle(G,F,E,blue, matte) );
	// triangles.push_back( Triangle(G,H,F,blue, matte) );


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

	vec3 centre;
	float radius;

	centre = vec3(0.45, 0.6, -0.1);
  radius = 0.4;
	spheres.push_back(Sphere(centre, radius, black, glass));

	centre = vec3(-0.5, 0.6, -0.8);
	radius = 0.3;
	spheres.push_back(Sphere(centre, radius, black, chrome));


}

#endif
