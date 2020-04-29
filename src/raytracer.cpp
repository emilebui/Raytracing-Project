// The JSON library allows you to reference JSON arrays like C++ vectors and JSON objects like C++ maps.

#include "raytracer.h"

#include <iostream>
#include <fstream>
#include <string>
#include <glm/glm.hpp>
#include <glm/gtx/string_cast.hpp>
#include <vector>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/matrix_inverse.hpp>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "json.hpp"

using json = nlohmann::json;

const char* PATH = "scenes/";

double fov = 60;
colour3 background_colour(0, 0, 0);
const float nAir = 1;
const int maxReflectLv = 8;
bool gothrough = false;
const double validT = 0.001;
const int PARTICLE_NUM = 20;
const int GLOSSY_INTERVAL = 2;
const int TEST_NUMS = 200;

//Test BVH
bool bvhScene = false;

json scene;

/****************************************************************************/

// here are some potentially useful utility functions

bool checkLightDis(point3 lightpos, point3 hitpt, point3 hitpt2);
bool checkHitPlane(point3 n, point3 a, point3 e, point3 d, point3& hitpt, float& t0);
float randF(float min, float max);
point3 getRandomPt(point3 p, float r);

//---------------------------------------------------------------
class Box {
public:
	float maxX;
	float minX;
	float maxY;
	float minY;
	float maxZ;
	float minZ;

	Box() {}

	Box(float x0, float x1, float y0, float y1, float z0, float z1)
	{
		minX = x0; maxX = x1; minY = y0; maxY = y1; minZ = z0; maxZ = z1;
	}

	bool checkIntersect(point3 d, point3 e) {
		float t0 = (minX - e.x) / d.x;
		float t1 = (maxX - e.x) / d.x;
		if (t0 > t1)
			std::swap(t0, t1);
		float tempT0 = (minY - e.y) / d.y;
		float tempT1 = (maxY - e.y) / d.y;
		if (tempT0 > tempT1)
			std::swap(tempT0, tempT1);
		if ((t0 > tempT1) || (tempT0 > t1))
			return false;
		if (tempT0 > t0)
			t0 = tempT0;
		if (tempT1 < t1)
			t1 = tempT1;
		tempT0 = (minZ - e.z) / d.z;
		tempT1 = (maxZ - e.z) / d.z;
		if (tempT0 > tempT1)
			std::swap(tempT0, tempT1);
		if ((t0 > tempT1) || (tempT0 > t1))
			return false;
		return true;
	}
};

class Object {
public:
	point3 amb;
	point3 diff;
	point3 spec;
	float shiny;
	point3 trans;
	point3 reflect;
	float refract;
	float oren;
	float glossy;
	bool toon = false;
	point3 c;
	Box b;

	virtual bool checkIntersect(point3 d, point3 e, point3& hitpt, point3& N, float& t0) = 0;
};

std::vector<Object*> objectList;

class Tri : public Object {

public:
	std::vector<point3> ps;

	void generateBox() {
		float minX, maxX, minY, maxY, minZ, maxZ;
		point3 tempP = ps[0];
		minX = tempP.x; maxX = tempP.x; minY = tempP.y; maxY = tempP.y; minZ = tempP.z; maxZ = tempP.z;

		for (int i = 0; i < ps.size(); i++) {
			tempP = ps[i];

			if (tempP.x < minX)
				minX = tempP.x;
			else if (tempP.x > maxX)
				maxX = tempP.x;

			if (tempP.y < minY)
				minY = tempP.y;
			else if (tempP.y > maxY)
				maxY = tempP.y;

			if (tempP.z < minZ)
				minZ = tempP.z;
			else if (tempP.z > maxZ)
				maxZ = tempP.z;
		}

		b = Box(minX, maxX, minY, maxY, minZ, maxZ);
	}

	Tri(std::vector<point3> points) {
		ps = points;

		float avgX = (ps[0].x + ps[1].x + ps[2].x) / 3;
		float avgY = (ps[0].y + ps[1].y + ps[2].y) / 3;
		float avgZ = (ps[0].z + ps[1].z + ps[2].z) / 3;
		c = point3(avgX, avgY, avgZ);
		generateBox();
	}
	void calculateNormal(point3& N) {
		point3 AB = ps[1] - ps[0];
		point3 AC = ps[2] - ps[0];
		N = glm::cross(AB, AC);
		N = glm::normalize(N);
	}

	void updatePts(point3 p, float s) {
		for (int i = 0; i < ps.size(); i++) {
			ps[i] = ps[i] * s;
			ps[i] = ps[i] + p;
		}
	}

	bool checkIntersect(point3 d, point3 e, point3& X, point3& N, float& t0) {
		calculateNormal(N);
		if (checkHitPlane(N, ps[0], e, d, X, t0)) {
			point3 BA = ps[1] - ps[0];
			point3 CB = ps[2] - ps[1];
			point3 AC = ps[0] - ps[2];
			point3 XA = X - ps[0];
			point3 XB = X - ps[1];
			point3 XC = X - ps[2];

			double tempA = glm::dot(glm::cross(BA, XA), N);
			double tempB = glm::dot(glm::cross(CB, XB), N);
			double tempC = glm::dot(glm::cross(AC, XC), N);

			if (tempA > 0 && tempB > 0 && tempC > 0)
			{
				return true;
			}
			else if (tempA < 0 && tempB < 0 && tempC < 0)
			{
				return true;
			}
		}
		return false;
	}
};

class Mesh : public Object {
public:
	std::vector<Tri> tris;
	void generateBox() {
		float minX, maxX, minY, maxY, minZ, maxZ;
		Tri TempT = tris[0];
		point3 center;
		float count = 0;
		std::vector<point3> pts = TempT.ps;
		minX = pts[0].x; maxX = pts[0].x;
		minY = pts[0].y; maxY = pts[0].y;
		minZ = pts[0].z; maxZ = pts[0].z+0.2;

		for (int i = 0; i < tris.size(); i++) {

			TempT = tris[i];
			pts = TempT.ps;

			for (int j = 0; j < pts.size(); j++) {
				point3 tempP = pts[j];

				if (tempP.x < minX)
					minX = tempP.x;
				if (tempP.x > maxX)
					maxX = tempP.x;
				if (tempP.y < minY)
					minY = tempP.y;
				if (tempP.y > maxY)
					maxY = tempP.y;
				if (tempP.z < minZ)
					minZ = tempP.z;
				if (tempP.z > maxZ)
					maxZ = tempP.z;
			}

			center += TempT.c;
			count++;
		}

		b = Box(minX, maxX, minY, maxY, minZ, maxZ);
		c = center / count;
	}

	Mesh(std::vector<Tri> tria) {
		tris = tria;
		generateBox();
	}
	bool checkIntersect(point3 d, point3 e, point3& hitpt, point3& N, float& t0) {
		for (int i = 0; i < tris.size(); i++) {
			Tri tempTri = tris[i];
			if (tempTri.checkIntersect(d, e, hitpt, N, t0)) {
				return true;
			}
		}
		return false;
	}
};

class Sphere : public Object {
public:
	float r;

	Box generateBox() {
		float minX, maxX, minY, maxY, minZ, maxZ;
		minX = c.x - r; maxX = c.x + r;
		minY = c.y - r; maxY = c.y + r;
		minZ = c.z - r; maxZ = c.z + r;
		return Box(minX, maxX, minY, maxY, minZ, maxZ);
	}

	Sphere(float rad, point3 p0) {
		r = rad;
		c = p0;
		b = generateBox();
	}

	bool checkIntersect(point3 d, point3 e, point3& hitpt, point3& N, float& t0) {
		double temp = pow(glm::dot(d, e - c), 2) - (glm::dot(d, d) * (glm::dot(e - c, e - c) - r * r));

		if (temp < 0)
			return false;

		double leftside = glm::dot(-d, e - c) / glm::dot(d, d);
		double t;

		if (temp > 0) {
			double rightside = sqrt(temp) / glm::dot(d, d);

			double t1 = leftside + rightside;
			double t2 = leftside - rightside;

			if (t1 < t2)
				t = t1;
			else
				t = t2;
		}
		else {
			t = leftside;
		}

		if (t < validT)
			return false;

		hitpt = e + float(t) * d;
		N = glm::normalize(hitpt - c);
		t0 = float(t);
		return true;
	}

};

class Plane : public Object {
public:
	point3 a;
	point3 n;

	Plane(point3 pos, point3 normal) {
		a = pos;
		n = normal;
	}
	bool checkIntersect(point3 d, point3 e, point3& hitpt, point3& N, float& t0) {
		N = glm::normalize(n);
		return checkHitPlane(n, a, e, d, hitpt, t0);
	}
};

std::vector<Tri> generateParticle(std::vector<Tri> tria, point3 pos, float r, float s) {
	std::vector<Tri> tris;

	for (int i = 0; i < PARTICLE_NUM; i++) {
		std::vector<Tri> temp;
		temp = tria;
		point3 newP = getRandomPt(pos, r);

		for (int j = 0; j < temp.size(); j++) {
			temp[j].updatePts(newP, s);
			tris.push_back(temp[j]);
		}
	}

	return tris;
}

class Particles : public Object {
public:
	std::vector<Tri> tris;

	Particles(std::vector<Tri> tria, point3 pos, float r, float s) {
		tris = generateParticle(tria, pos, r, s);
	}

	bool checkIntersect(point3 d, point3 e, point3& hitpt, point3& N, float& t0) {
		for (int i = 0; i < tris.size(); i++) {
			Tri tempTri = tris[i];
			if (tempTri.checkIntersect(d, e, hitpt, N, t0)) {
				return true;
			}
		}
		return false;
	}


};

class Light {
public:
	point3 color;
	bool amb;
	virtual point3 getL(point3 hitpt) { return point3(0); }
	virtual bool checkLightDist(point3 hitpt, point3 hitpt2) { return false; }
	virtual bool checkValid(point3 L1) { return false; }
};

class PointLight : public Light {
public:
	point3 p;
	PointLight(point3 pos, point3 c) {
		p = pos;
		color = c;
	}
	bool checkLightDist(point3 hitpt, point3 hitpt2) {
		return checkLightDis(p, hitpt, hitpt2);
	}
	point3 getL(point3 hitpt) {
		point3 L = glm::normalize(p - hitpt);
		return L;
	}
	bool checkValid(point3 L1) {
		return true;
	}
};

class SpotLight : public Light {
public:
	point3 p;
	point3 dir;
	float angle;
	SpotLight(point3 pos, point3 d, float cutoff, point3 c) {
		p = pos;
		dir = d;
		angle = cutoff;
		color = c;
	}
	bool checkLightDist(point3 hitpt, point3 hitpt2) {
		return checkLightDis(p, hitpt, hitpt2);
	}
	point3 getL(point3 hitpt) {
		point3 L = glm::normalize(p - hitpt);
		return L;
	}
	bool checkValid(point3 L1) {
		point3 D = glm::normalize(dir);
		float tempCosine = glm::dot(D, L1);
		if (tempCosine >= cos(glm::radians(angle)))
			return true;

		return false;
	}
};

class DirectionalLight : public Light {
public:
	point3 dir;
	DirectionalLight(point3 d, point3 c) {
		dir = d;
		color = c;
	}
	bool checkLightDist(point3 hitpt, point3 hitpt2) {
		return true;
	}
	point3 getL(point3 hitpt) {
		point3 L = glm::normalize(-dir);
		return L;
	}
	bool checkValid(point3 L1) {
		return true;
	}
};

class AmbientLight : public Light {
public:
	AmbientLight(point3 c) {
		color = c;
		amb = true;
	}
};
//--------------------------------------------------------------
//BVH
//--------------------------------------------------------------
class BVHNode {
public:
	BVHNode* leftnode;
	BVHNode* rightnode;
	Object* leftobj = NULL;
	Object* rightobj = NULL;
	Box b;
	bool bottom = false;


	void makeBox(std::vector<Object*> list) {
		float minX, maxX, minY, maxY, minZ, maxZ;

		//setup
		Box tempB = list[0]->b;
		minX = tempB.minX; maxX = tempB.maxX; minY = tempB.minY; maxY = tempB.maxY; minZ = tempB.minZ; maxZ = tempB.maxZ;

		for (int i = 0; i < list.size(); i++) {
			tempB = list[i]->b;

			if (tempB.minX < minX)
				minX = tempB.minX;
			if (tempB.maxX > maxX)
				maxX = tempB.maxX;
			if (tempB.minY < minY)
				minY = tempB.minY;
			if (tempB.maxY > maxY)
				maxY = tempB.maxY;
			if (tempB.minZ < minZ)
				minZ = tempB.minZ;
			if (tempB.maxZ > maxZ)
				maxZ = tempB.maxZ;
		}

		b = Box(minX, maxX, minY, maxY, minZ, maxZ);
	}

	void makeBVHNode(std::vector<Object*> list) {
		
		//setup
		std::vector<Object*> leftL;
		std::vector<Object*> rightL;

		//Choose axis
		float lengthX = b.maxX - b.minX;
		float lengthY = b.maxY - b.minY;
		float lengthZ = b.maxZ - b.minZ;
		float midP;

		if (lengthX >= lengthY && lengthX >= lengthZ) {

			midP = (b.maxX + b.minX) / 2.0f;

			for (int i = 0; i < list.size(); i++) {
				point3 center = list[i]->c;

				if (center.x >= midP) {
					rightL.push_back(list[i]);
				}
				else {
					leftL.push_back(list[i]);
				}
			}
		}
		else if (lengthY >= lengthX && lengthY >= lengthZ) {
			midP = (b.maxY + b.minY) / 2.0f;

			for (int i = 0; i < list.size(); i++) {
				point3 center = list[i]->c;

				if (center.y >= midP) {
					rightL.push_back(list[i]);
				}
				else {
					leftL.push_back(list[i]);
				}
			}
		}
		else if (lengthZ >= lengthX && lengthZ >= lengthY) {
			midP = (b.maxZ + b.minZ) / 2.0f;

			for (int i = 0; i < list.size(); i++) {
				point3 center = list[i]->c;

				if (center.z >= midP) {
					rightL.push_back(list[i]);
				}
				else {
					leftL.push_back(list[i]);
				}
			}
		}

		leftnode = new BVHNode(leftL);
		rightnode = new BVHNode(rightL);
	}

	BVHNode(std::vector<Object*> list) {
		
		if (list.size() <= 2) {
			bottom = true;

			if (list.size() == 2) {
				leftobj = list[0];
				rightobj = list[1];
			}
			else if (list.size() == 1)
				leftobj = list[0];
		}
		else {
			makeBox(list);
			makeBVHNode(list);
		}
	}

	Object* checkIntersect(point3 d, point3 e, point3& hitpt, point3& N, float& t0, bool& check) {

		if (!bottom) {
			check = b.checkIntersect(d,e);
			
			if (check) {
				bool checkL;
				bool checkR;
				float tempTL;
				float tempTR;
				Object* tempOL = leftnode->checkIntersect(d,e,hitpt,N, tempTL,checkL);
				Object* tempOR = rightnode->checkIntersect(d, e, hitpt, N, tempTR, checkR);

				if (checkL) {
					if (checkR) {
						check = true;
						if (tempTL <= tempTR) {
							t0 = tempTL;
							return tempOL;
						}
						else {
							t0 = tempTR;
							return tempOR;
						}
							
					}
					else {
						t0 = tempTL;
						check = true;
						return tempOL;
					}
				}
				else if (checkR) {
					t0 = tempTR;
					check = true;
					return tempOR;
				}
				else {
					check = false;
					return NULL;
				}

			}
			else {
				return NULL;
			}
		}
		else {
			float tempTL;
			float tempTR;

			if (leftobj != NULL && leftobj->checkIntersect(d, e, hitpt, N, tempTL)) {
				if (rightobj != NULL && rightobj->checkIntersect(d, e, hitpt, N, tempTR)) {
					check = true;
					if (tempTL <= tempTR) {
						t0 = tempTL;
						return leftobj;
					}
					else {
						t0 = tempTR;
						return rightobj;
					}
						
				}
				else {
					t0 = tempTL;
					check = true;
					return leftobj;
				}
			}
			else if (rightobj != NULL && rightobj->checkIntersect(d, e, hitpt, N, tempTR)) {
				t0 = tempTR;
				check = true;
				return rightobj;
			}
			else {
				check = false;
				return NULL;
			}
		}
		
	}
};

class BVH {
public:
	BVHNode* leftnode;
	BVHNode* rightnode;

	BVH() {}

	BVH(std::vector<Object*> list) {
		std::vector<Object*> leftL;
		std::vector<Object*> rightL;

		for (int i = 0; i < list.size(); i++) {
			point3 center = list[i]->c;

			if (center.x >= 0) {
				rightL.push_back(list[i]);
			}
			else {
				leftL.push_back(list[i]);
			}
		}

		leftnode = new BVHNode(leftL);
		rightnode = new BVHNode(rightL);
	}
};

Object* traceObjectBVH(BVH bv, point3 d, point3 e, point3& hitpt, point3& N, float& t, bool& check) {

	bool checkL;
	bool checkR;
	float tempTL;
	float tempTR;
	Object* tempOL = bv.leftnode->checkIntersect(d, e, hitpt, N, tempTL, checkL);
	Object* tempOR = bv.rightnode->checkIntersect(d, e, hitpt, N, tempTR, checkR);

	if (checkL) {
		if (checkR) {
			check = true;
			if (tempTL <= tempTR) {
				t = tempTL;
				return tempOL;
			}
			else {
				t = tempTR;
				return tempOR;
			}

		}
		else {
			check = true;
			t = tempTL;
			return tempOL;
		}
	}
	else if (checkR) {
		check = true;
		t = tempTR;
		return tempOR;
	}
	else {
		check = false;
		return NULL;
	}
	

}


//--------------------------------------------------------------
//Prototype
//--------------------------------------------------------------
bool rayTrace(const point3& e, const point3& s, point3& colour, bool doColor, int level);
bool checkIfInShadow(point3 L2, point3 hitpt, int level, Light* light);
glm::vec3 vector_to_vec3(const std::vector<float>& v);
Object* traceObjectAll(point3 d, point3 e, bool& check, point3& hitpt0);

std::vector<Light*> lightList;
std::vector<Object*> planes;
BVH bvhList;
//---------------------------------------------------------------------------------------
std::vector<Tri> getTriList(std::vector<std::vector<std::vector<float>>> tris) {
	std::vector<Tri> triList;
	for (int i = 0; i < tris.size(); i++) {
		std::vector<point3> ps;
		std::vector<std::vector<float>> tri = tris[i];
		point3 A = vector_to_vec3(tri[0]);
		point3 B = vector_to_vec3(tri[1]);
		point3 C = vector_to_vec3(tri[2]);
		ps.push_back(A);
		ps.push_back(B);
		ps.push_back(C);
		triList.push_back(Tri(ps));
	}

	return triList;
}

float randF(float min, float max) {
	//srand(time(NULL));
	float percentage = ( (float)rand() ) / (float)RAND_MAX;
	float d = max - min;
	float r = percentage * d;
	float result = min + r;
	return result;
}

point3 getRandomPt(point3 p, float r) {
	float x = randF(p.x - r, p.x + r);
	float y = randF(p.y - r, p.y + r);
	float z = randF(p.z - r, p.z + r);
	return point3(x, y, z);
}

void copyStats(Object* tempO, Object* o) {
	tempO->amb = o->amb;
	tempO->diff = o->diff;
	tempO->glossy = o->glossy;
	tempO->oren = o->oren;
	tempO->reflect = o->reflect;
	tempO->refract = o->refract;
	tempO->shiny = o->shiny;
	tempO->spec = o->spec;
	tempO->toon = o->toon;
	tempO->trans = o->trans;
}

//---------------------------------------------------------------------------------------
void makeObjectList() {
	json objects = scene["objects"];

	for (json::iterator it = objects.begin(); it != objects.end(); ++it) {
		json& object = *it;
		Object* o;
		bool triangles = false;
		bool plane = false;
		std::vector<Tri> triList;
		if (object["type"] == "sphere") {
			std::vector<float> pos = object["position"];
			point3 c = vector_to_vec3(pos);
			float r = float(object["radius"]);
			o = new Sphere(r, c);
		}
		else if (object["type"] == "plane") {
			std::vector<float> pos = object["position"];
			point3 a = vector_to_vec3(pos);
			std::vector<float> normal = object["normal"];
			point3 n = vector_to_vec3(normal);
			o = new Plane(a, n);
			plane = true;
		}
		else if (object["type"] == "mesh") {
			std::vector<std::vector<std::vector<float>>> tris = object["triangles"];
			triList = getTriList(tris);
			o = new Mesh(triList);
			//triangles = true;
		}
		else if (object["type"] == "particles") {
			std::vector<std::vector<std::vector<float>>> tris = object["triangles"];
			triList = getTriList(tris);
			std::vector<float> pos = object["position"];
			point3 a = vector_to_vec3(pos);
			float r = float(object["radius"]);
			float s = float(object["scale"]);
			o = new Particles(triList, a, r, s);
			triangles = true;
			triList = generateParticle(triList, a, r, s);
		}

		//Get color
		json& material = object["material"];
		json test = material["ambient"];
		if (test.is_array()) {
			std::vector<float> ambient = material["ambient"];
			point3 ambRate = vector_to_vec3(ambient);
			o->amb = ambRate;
		}

		test = material["diffuse"];
		if (test.is_array()) {
			std::vector<float> tempV = material["diffuse"];
			point3 tempP = vector_to_vec3(tempV);
			o->diff = tempP;
		}

		test = material["specular"];
		if (test.is_array()) {
			std::vector<float> tempV = material["specular"];
			point3 tempP = vector_to_vec3(tempV);
			o->spec = tempP;
			float sh = material["shininess"];
			o->shiny = sh;
		}

		test = material["reflective"];
		if (test.is_array()) {
			std::vector<float> tempV = material["reflective"];
			point3 tempP = vector_to_vec3(tempV);
			o->reflect = tempP;
		}

		test = material["transmissive"];
		if (test.is_array()) {
			std::vector<float> tempV = material["transmissive"];
			point3 tempP = vector_to_vec3(tempV);
			o->trans = tempP;
		}
		test = material["refraction"];
		if (!test.is_null()) {
			float rf = material["refraction"];
			o->refract = rf;
		}
		test = material["oren"];
		if (!test.is_null()) {
			float ore = material["oren"];
			o->oren = ore;
		}
		test = material["glossy"];
		if (!test.is_null()) {
			float gl = material["glossy"];
			o->glossy = gl;
		}
		test = material["toon"];
		if (!test.is_null()) {
			bool to = material["toon"];
			o->toon = to;
		}

		//add
		if (!triangles && !plane)
			objectList.push_back(o);
		else if (triangles) {
			for (int i = 0; i < triList.size(); i++) {
				Tri tempTri = triList[i];
				Object* tempO = new Tri(tempTri.ps);
				copyStats(tempO, o);
				objectList.push_back(tempO);
			}
		}
		else if (plane) {
			planes.push_back(o);
		}
	}
}

void makeLightList() {
	json& lights = scene["lights"];

	for (json::iterator it = lights.begin(); it != lights.end(); ++it) {
		json& light = *it;
		Light* l;

		std::vector<float> tempC = light["color"];
		point3 color = vector_to_vec3(tempC);

		if (light["type"] == "ambient") {
			l = new AmbientLight(color);
		}
		else if (light["type"] == "point") {
			std::vector<float> tempV = light["position"];
			point3 tempP = vector_to_vec3(tempV);
			l = new PointLight(tempP, color);
		}
		else if (light["type"] == "directional") {
			std::vector<float> tempV = light["direction"];
			point3 tempP = vector_to_vec3(tempV);
			l = new DirectionalLight(tempP, color);
		}
		else if (light["type"] == "spot") {
			std::vector<float> tempV = light["position"];
			point3 pos = vector_to_vec3(tempV);
			std::vector<float> tempV2 = light["direction"];
			point3 dir = vector_to_vec3(tempV2);
			float cut = light["cutoff"];
			l = new SpotLight(pos, dir, cut, color);
		}

		lightList.push_back(l);
	}
}



//---------------------------------------------------------------------------------------
json find(json& j, const std::string key, const std::string value) {
	json::iterator it;
	for (it = j.begin(); it != j.end(); ++it) {
		if (it->find(key) != it->end()) {
			if ((*it)[key] == value) {
				return *it;
			}
		}
	}
	return json();
}

glm::vec3 vector_to_vec3(const std::vector<float>& v) {
	return glm::vec3(v[0], v[1], v[2]);
}

/****************************************************************************/

void choose_scene(char const* fn) {
	if (fn == NULL) {
		std::cout << "Using default input file " << PATH << "p.json\n";
		fn = "bvh";
	}

	std::cout << "Loading scene " << fn << std::endl;

	std::string fname = PATH + std::string(fn) + ".json";
	std::fstream in(fname);
	if (!in.is_open()) {
		std::cout << "Unable to open scene file " << fname << std::endl;
		exit(EXIT_FAILURE);
	}

	in >> scene;

	if (fn == "bvh")
		bvhScene = true;

	makeObjectList();
	makeLightList();

	if (bvhScene) {
		Object* chosenObj = objectList[0];
		//printf("%f \n", chosenObj->diff.x);


		for (int i = 0; i < TEST_NUMS; i++) {
			point3 randomP = getRandomPt(chosenObj->c, 2.0f);
			//printf("%f , %f , %f \n", randomP.x, randomP.y, randomP.z);
			Object* newS = new Sphere(0.10, randomP);
			newS->amb = chosenObj->amb;
			newS->diff = chosenObj->diff;
			newS->spec = chosenObj->spec;
			newS->shiny = chosenObj->shiny;
			newS->c = randomP;
			objectList.push_back(newS);
		}
	}

	bvhList = BVH(objectList);

	//printf("%d \n", planes.size());

	json camera = scene["camera"];
	// these are optional parameters (otherwise they default to the values initialized earlier)
	if (camera.find("field") != camera.end()) {
		fov = camera["field"];
		std::cout << "Setting fov to " << fov << " degrees.\n";
	}
	if (camera.find("background") != camera.end()) {
		background_colour = vector_to_vec3(camera["background"]);
		std::cout << "Setting background colour to " << glm::to_string(background_colour) << std::endl;
	}
}

//------------------------------------------------------------------------------------------------

bool checkHitPlane(point3 n, point3 a, point3 e, point3 d, point3& hitpt, float& t0) {

	point3 ae = a - e;
	double t = glm::dot(n, ae) / glm::dot(n, d);

	if (t > validT)
	{
		hitpt = e + float(t) * d;
		t0 = float(t);
		return true;
	}

	return false;
}
//-----------------------------------------------------------------------------------------
bool checkLightDis(point3 lightpos, point3 hitpt, point3 hitpt2) {

	float disL = glm::distance(lightpos, hitpt);
	float disO = glm::distance(hitpt2, hitpt);

	if (disL <= disO)
		return false;

	return true;
}



bool tracebackForShadow(point3 L2, point3 hitpt, Light* light) {
	bool check;
	float tempT;
	point3 tempP;
	point3 tempN;

	/**
	for (int i = 0; i < objectList.size(); i++) {
		Object* o = objectList[i];
		check = o->checkIntersect(L2, hitpt, tempP, tempN, tempT);
		if (check && light->checkLightDist(hitpt, tempP))
			return true;
	}
	**/
	
	point3 s2 = hitpt + L2;
	traceObjectAll(L2, hitpt, check, tempP);
	if (check && light->checkLightDist(hitpt, tempP))
		return true;


	return false;
}


bool checkIfInShadow(point3 L2, point3 hitpt, Light* light) {

	point3 tempC;

	bool check = tracebackForShadow(L2, hitpt, light);
	//bool check = false;

	return check;
}

point3 getTotalColorR(point3 c, point3 hitpt, int level, float glossy) {
	float x = c.x;
	float y = c.y;
	float z = c.z;
	point3 tempC;
	point3 totalC;
	point3 tempS;

	float invG = glossy / GLOSSY_INTERVAL;

	float count = 0;

	//+
	for (int i = 0; i < GLOSSY_INTERVAL; i++)
	{
		y = c.y;
		z = c.z;
		for (int j = 0; j < GLOSSY_INTERVAL; j++) 
		{
			
			z = c.z;
			for (int k = 0; k < GLOSSY_INTERVAL; k++) 
			{
				tempS = point3(x, y, z);
				if (!rayTrace(hitpt, tempS, tempC, true, level)) {
					tempC = background_colour;
				}
				totalC += tempC;
				count++;
				z += invG;
			}
			y += invG;
		}
		x += invG;
	}

	//-
	x = c.x;
	for (int i = 0; i < GLOSSY_INTERVAL; i++)
	{
		y = c.y;
		z = c.z;
		for (int j = 0; j < GLOSSY_INTERVAL; j++)
		{
			z = c.z;
			for (int k = 0; k < GLOSSY_INTERVAL; k++)
			{
				tempS = point3(x, y, z);
				if (!rayTrace(hitpt, tempS, tempC, true, level)) {
					tempC = background_colour;
				}
				totalC += tempC;
				count++;
				z -= invG;
			}
			y -= invG;
		}
		x -= invG;
	}

	point3 avgC = totalC / count;

	return avgC;
}



point3 getReflectColor(point3 R, point3 hitpt, int level, float glossy) {

	point3 resultColor = point3(0);

	point3 R2 = hitpt + R;

	if (glossy > 0) {
		resultColor = getTotalColorR(R2, hitpt, level, glossy);
	}
	else {
		if (!rayTrace(hitpt, R2, resultColor, true, level))
		{
			resultColor = background_colour;
		}
	}

	

	return resultColor;

}

point3 getTransparentColor(point3 R, point3 hitpt, int level) {
	point3 resultColor = point3(0);
	point3 R2 = hitpt + R;

	//tracebackForShadow(R2, hitpt, true, resultColor, level);
	if (!rayTrace(hitpt, R2, resultColor, true, level))
	{
		resultColor = background_colour;
	}


	return resultColor;
}
//-----------------------------------------------------------------------------------------
float calculateOrenColor(point3 N, point3 L, point3 V, float oren) {

	float oren2 = pow(oren, 2);

	float NL = glm::dot(N, L);
	float NV = glm::dot(N, V);
	float Oi = acos(NL);
	float Or = acos(NV);

	float tempC = std::max(glm::dot(glm::normalize(V - N * NV), glm::normalize(L - N * NL) ), 0.0f );
	float alp = std::max(Oi, Or);
	float bet = std::min(Oi, Or);
	float A = 1.0f - (0.5f * oren2 / (oren2 + 0.33f));
	float B = 0.45f * oren2 / (oren2 + 0.09);
	float tempA = std::max(glm::dot(N, L), 0.0f);

	float result = tempA * (A + B * tempC * sin(alp) * tan(bet));
	return result;
}

//------------------------------------------------------------------------------------------
point3 calculateLightColor(Light* light, point3& L, point3& N, point3& V, Object* object, point3 hitpt) {
	point3 result = point3(0);

	point3 Ids = light->color;
	bool toon_shading = object->toon;

	point3 diffRate = object->diff;
	if (diffRate.x > 0 || diffRate.y > 0 || diffRate.z > 0) {
		if (object->oren > 0) {
			float tempA = calculateOrenColor(N, L, V, object->oren);
			result = tempA * Ids * diffRate;
		}
		else {
			if (toon_shading) {
				float tempA = glm::dot(N, L);
				/**if (tempA > 0) {
					tempA = 1;
				}
				else {
					tempA = 0;
				}**/
				tempA = glm::smoothstep(0.0f, 0.01f, tempA);
				result = Ids * diffRate * tempA;
			}
			else {
				float tempA = std::max(glm::dot(N, L), 0.0f);
				result = Ids * diffRate * tempA;
			}
			
		}
		
	}


	point3 specular = object->spec;
	float shiny = object->shiny;

	point3 R;

	if (shiny > 0.01) {
		R = 2 * glm::dot(N, L) * N - L;
		float angle = std::max(glm::dot(R, V), 0.0f);
		float spec = pow(angle, shiny);

		if (toon_shading) {
			spec = glm::smoothstep(0.005f, 0.01f, spec);
		}

		//float angle = glm::dot(R, V);

		point3 tempR = Ids * specular * spec;

		result += tempR;
	}

	return result;
}

//------------------------------------------------------------------------------------------
point3 calculateLight(Light* light, point3& hitpt, point3& N, point3& V, Object* object) {

	point3 result = point3(0);

	if (light->amb == true) {
		result = light->color * object->amb;
	}
	else {
		point3 L = light->getL(hitpt);

		if (!checkIfInShadow(L, hitpt, light) && light->checkValid(-L))
		{
			result = calculateLightColor(light, L, N, V, object, hitpt);
		}

	}
	return result;
}

//------------------------------------------------------------------------------------------
void calculateColor(Object* object, point3& hitpt, point3& N, point3& V, point3& colour, int level) {

	point3 tempC;
	point3 tempC2;
	point3 trueColor;
	point3 reflectionColor;
	point3 transparentColor;

	for (int i = 0; i < lightList.size(); i++) {
		Light* light = lightList[i];
		tempC2 = calculateLight(light, hitpt, N, V, object);
		tempC2 = glm::clamp(tempC2, 0.0f, 1.0f);
		tempC += tempC2;
	}

	trueColor = glm::clamp(tempC, 0.0f, 1.0f);

	point3 reflective = object->reflect;
	if (reflective.x > 0.01 && level < maxReflectLv)
	{
		//float angle = std::max(glm::dot(N, V), 0.0f);
		point3 R = 2 * glm::dot(N, V) * N - V;
		//point3 R = 2 * angle * N - V;
		R = glm::normalize(R);
		reflectionColor = getReflectColor(R, hitpt, level + 1, object->glossy);
		reflectionColor = reflectionColor * reflective;
		reflectionColor = glm::clamp(reflectionColor, 0.0f, 1.0f);
		//tempC += newColor;
	}

	transparentColor = trueColor;

	point3 transparent = object->trans;
	if (transparent.x > 0.01) {
		bool refract = false;
		point3 V2 = -V;
		V2 = glm::normalize(V2);

		float refractRate = object->refract;
		if (refractRate > 0.01)
		{
			point3 n = glm::normalize(N);
			float angle = glm::clamp(glm::dot(V2, n), -1.0f, 1.0f);
			//float angle = glm::dot(V2, n);
			float nI;
			float nR;
			gothrough = !gothrough;

			if (angle < 0) {

				nI = nAir;
				nR = refractRate;
				angle = -angle;
				//n = -n;
			}
			else {
				nI = refractRate;
				nR = nAir;
				n = -n;
			}

			//float tempF = 1 - (((nI * nI) * (1 - pow(angle, 2))) / (nR * nR));
			float tempF = 1 - nR * nR * (1 - angle * angle);
			point3 vR;
			//pow(angle, 2) <= (1 - (nI * nI / (nR * nR)))
			if (tempF >= 0)
			{
				//vR = ((nI * (V2 - (n * angle))) / nR) - (n * sqrt(tempF));
				float nIR = nI / nR;
				vR = nIR * V2 + (nIR * angle - sqrt(tempF)) * n;
				V2 = vR;
				V2 = glm::normalize(V2);
			}
			else {
				refract = true;
			}
		}

		if (!refract) {
			point3 nColor = getTransparentColor(V2, hitpt, level);
			nColor = nColor * transparent;
			nColor = glm::clamp(nColor, 0.0f, 1.0f);
			point3 one = point3(1, 1, 1);

			transparentColor = transparentColor * (one - transparent) + nColor;
		}
		else {
			point3 R2 = 2 * glm::dot(N, V) * N - V;
			point3 nColor = getReflectColor(R2, hitpt, level + 1, object->glossy);
			nColor = nColor * transparent;
			nColor = glm::clamp(nColor, 0.0f, 1.0f);
			tempC += nColor;
		}


		//tempC += nColor;
	}

	tempC = transparentColor + reflectionColor;
	colour = glm::clamp(tempC, 0.0f, 1.0f);

	//colour = tempC;
}

//------------------------------------------------------------------------------------------
Object* traceObjectAll(point3 d, point3 e, bool& check, point3& hitpt0) {
	Object* chosenObj;
	point3 hitpt;
	point3 N;
	bool check2;
	bool checkB;
	float tempT1;
	float tempT;
	float minT = 99;

	chosenObj = traceObjectBVH(bvhList, d, e, hitpt, N, minT, check);
	if (!check) {
		minT = 99;
	}
	hitpt0 = hitpt;
	//Check planes
	for (int i = 0; i < planes.size(); i++)
	{
		Object* o = planes[i];
		check2 = o->checkIntersect(d, e, hitpt, N, tempT);
		if (check2 && tempT < minT)
		{
			minT = tempT;
			chosenObj = o;
			check = true;
			hitpt0 = hitpt;
		}
	}

	return chosenObj;
}


bool rayTrace(const point3& e, const point3& s, point3& colour, bool doColor, int level) {

	bool check = false;
	float tempT;
	float minT = 99;
	Object* chosenObj;
	point3 hitpt;
	point3 N;
	bool check2;
	point3 d = s - e;

	chosenObj = traceObjectAll(d, e, check, hitpt);
	

	if (check) {
		if (doColor)
		{
			chosenObj->checkIntersect(d, e, hitpt, N, tempT);
			point3 V = glm::normalize(e - hitpt);
			calculateColor(chosenObj, hitpt, N, V, colour, level);
		}
		return true;
	}

	return false;
}




bool trace(const point3& e, const point3& s, colour3& colour, bool pick) {
	// NOTE 1: This is a demo, not ray tracing code! You will need to replace all of this with your own code...
  // NOTE 2: You can work with JSON objects directly (like this sample code), read the JSON objects into your own data structures once and render from those (probably in choose_scene), or hard-code the objects in your own data structures and choose them by name in choose_scene; e.g. choose_scene('e') would pick the same scene as the one in "e.json". Your choice.
  // If you want to use this JSON library, https://github.com/nlohmann/json for more information. The code below gives examples of everything you should need: getting named values, iterating over arrays, and converting types.

	// traverse the objects
	return rayTrace(e, s, colour, true, 0);
}