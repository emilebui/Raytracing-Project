// Assignment 2 Question 1

#include "common.h"
#include "raytracer.h"

#include <iostream>
#define M_PI 3.14159265358979323846264338327950288
#include <cmath>

#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/matrix_inverse.hpp>
#include <glm/gtc/type_ptr.hpp>

const char *WINDOW_TITLE = "Ray Tracing";
const double FRAME_RATE_MS = 1;

colour3 texture[1<<16]; // big enough for a row of pixels
point3 vertices[2]; // xy+u for start and end of line
GLuint Window;
int vp_width, vp_height;
float drawing_y = 0;


bool antialiasing = false;

point3 eye;
float d = 1;

bool cam = false;
//CAMERA POSITION (CHOOSE from 0 - 3)
const int campos = 0;


//camera vars
float angle;
glm::vec3 rotateCoord;
point3 newpos;

//----------------------------------------------------------------------------

point3 s(int x, int y) {
	float aspect_ratio = (float)vp_width / vp_height;
	float h = d * (float)tan((M_PI * fov) / 180.0 / 2.0);
	float w = h * aspect_ratio;
   
	float top = h;
	float bottom = -h;
	float left = -w;
	float right = w;
   
	float u = left + (right - left) * (x + 0.5f) / vp_width;
	float v = bottom + (top - bottom) * (y + 0.5f) / vp_height;
   
	return point3(u, v, -d);
}
//
void chooseCameraPosition() {
	if (campos == 0) {
		angle = glm::radians(30.0f);
		rotateCoord = glm::vec3(0, 1, 0);
		newpos = point3(2, 0, 0);
	}
	else if (campos == 1) {
		angle = glm::radians(20.0f);
		rotateCoord = glm::vec3(0, 0, 1);
		newpos = point3(0, 2, 0);
	}
	else if (campos == 2) {
		angle = glm::radians(50.0f);
		rotateCoord = glm::vec3(1, 0, 1);
		newpos = point3(0, 0, 1);
	}
	else {
		angle = glm::radians(90.0f);
		rotateCoord = glm::vec3(1, 0, 0);
		newpos = point3(0, 1, 1);
	}


}

//----------------------------------------------------------------------------

// OpenGL initialization
void init(char *fn) {
	choose_scene(fn);
   
	// Create a vertex array object
	GLuint vao;
	glGenVertexArrays( 1, &vao );
	glBindVertexArray( vao );

	// Create and initialize a buffer object
	GLuint buffer;
	glGenBuffers( 1, &buffer );
	glBindBuffer( GL_ARRAY_BUFFER, buffer );
	glBufferData( GL_ARRAY_BUFFER, sizeof(vertices), NULL, GL_STATIC_DRAW );

	// Load shaders and use the resulting shader program
	GLuint program = InitShader( "v.glsl", "f.glsl" );
	glUseProgram( program );

	// set up vertex arrays
	GLuint vPos = glGetAttribLocation( program, "vPos" );
	glEnableVertexAttribArray( vPos );
	glVertexAttribPointer( vPos, 3, GL_FLOAT, GL_FALSE, 0, 0 );

	Window = glGetUniformLocation( program, "Window" );

	// glClearColor( background_colour[0], background_colour[1], background_colour[2], 1 );
	glClearColor( 0.7, 0.7, 0.8, 1 );

	// set up a 1D texture for each scanline of output
	GLuint textureID;
	glGenTextures( 1, &textureID );
	glBindTexture( GL_TEXTURE_1D, textureID );
	glTexParameteri( GL_TEXTURE_1D, GL_TEXTURE_MAG_FILTER, GL_NEAREST );
	glTexParameteri( GL_TEXTURE_1D, GL_TEXTURE_MIN_FILTER, GL_NEAREST );
}

//----------------------------------------------------------------------------

void display( void ) {
	// draw one scanline at a time, to each buffer; only clear when we draw the first scanline
	// (when fract(drawing_y) == 0.0, draw one buffer, when it is 0.5 draw the other)
	
	if (drawing_y <= 0.5) {
		glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );

		glFlush();
		glFinish();
		glutSwapBuffers();

		drawing_y += 0.5;

	} else if (drawing_y >= 1.0 && drawing_y <= vp_height + 0.5) {
		int y = int(drawing_y) - 1;

		// only recalculate if this is a new scanline
		if (drawing_y == int(drawing_y)) {

			for (int x = 0; x < vp_width; x++) {

				point3 vs = s(x, y);

				//Rotate and translate cam
				if (cam) {

					chooseCameraPosition();

					//Rotate + transition camera
					glm::mat4 matrix, trans, rot;
					rot = glm::rotate(rot, angle, rotateCoord);
					trans = glm::translate(glm::mat4(), newpos);
					matrix = trans * rot;

					glm::vec4 vs2 = glm::vec4(vs.x, vs.y, vs.z, 0);
					vs2 = vs2 * matrix;
					vs = point3(vs2.x, vs2.y, vs2.z);
					glm::vec4 eye2 = glm::vec4(eye.x, eye.y, eye.z, 0);
					eye2 = eye2 * matrix;
					eye = point3(eye2.x, eye2.y, eye2.z);
				}
				
				//Antialiasing
				if (!antialiasing) {
					if (!trace(eye, vs, texture[x], false)) {
						texture[x] = background_colour;
					}
				}
				else {
					point3 tempC;
					point3 totalC;
					point3 tempS;

					//first ray
					if (!trace(eye, vs, tempC, false)) {
						tempC = background_colour;
					}
					totalC += tempC;

					//second ray
					tempS = point3(vs.x, vs.y + 0.001, vs.z);
					if (!trace(eye, tempS, tempC, false)) {
						tempC = background_colour;
					}
					totalC += tempC;

					//third ray
					tempS = point3(vs.x + 0.001, vs.y, vs.z);
					if (!trace(eye, tempS, tempC, false)) {
						tempC = background_colour;
					}
					totalC += tempC;

					//forth ray
					tempS = point3(vs.x - 0.001, vs.y, vs.z);
					if (!trace(eye, tempS, tempC, false)) {
						tempC = background_colour;
					}
					totalC += tempC;

					//fifth ray
					tempS = point3(vs.x, vs.y - 0.001, vs.z);
					if (!trace(eye, tempS, tempC, false)) {
						tempC = background_colour;
					}
					totalC += tempC;

					point3 avgC = totalC / 5.0f;
					texture[x] = avgC;

				}



				
			}

			// to ensure a power-of-two texture, get the next highest power of two
			// https://graphics.stanford.edu/~seander/bithacks.html#RoundUpPowerOf2
			unsigned int v; // compute the next highest power of 2 of 32-bit v
			v = vp_width;
			v--;
			v |= v >> 1;
			v |= v >> 2;
			v |= v >> 4;
			v |= v >> 8;
			v |= v >> 16;
			v++;
			
			glTexImage1D( GL_TEXTURE_1D, 0, GL_RGB, v, 0, GL_RGB, GL_FLOAT, texture );
			vertices[0] = point3(0, y, 0);
			vertices[1] = point3(v, y, 1);
			glBufferSubData( GL_ARRAY_BUFFER, 0, 2 * sizeof(point3), vertices);
		}

		glDrawArrays( GL_LINES, 0, 2 );
		
		glFlush();
		glFinish();
		glutSwapBuffers();
		
		drawing_y += 0.5;
	}
}

//----------------------------------------------------------------------------

void keyboard( unsigned char key, int x, int y ) {
	switch( key ) {
	case 033: // Escape Key
	case 'q': case 'Q':
		exit( EXIT_SUCCESS );
		break;
	case ' ':
		drawing_y = 1;
		break;
	}
}

//----------------------------------------------------------------------------

void mouse( int button, int state, int x, int y ) {
	y = vp_height - y - 1;
	if ( state == GLUT_DOWN ) {
		switch( button ) {
		case GLUT_LEFT_BUTTON:
			colour3 c;
			point3 uvw = s(x, y);
			std::cout << std::endl;
			if (trace(eye, uvw, c, true)) {
				std::cout << "HIT @ ( " << uvw.x << "," << uvw.y << "," << uvw.z << " )\n";
				std::cout << "      colour = ( " << c.r << "," << c.g << "," << c.b << " )\n";
			} else {
				std::cout << "MISS @ ( " << uvw.x << "," << uvw.y << "," << uvw.z << " )\n";
			}
			break;
		}
	}
}

//----------------------------------------------------------------------------

void update( void ) {
}

//----------------------------------------------------------------------------

void reshape( int width, int height ) {
	glViewport( 0, 0, width, height );

	// GLfloat aspect = GLfloat(width)/height;
	// glm::mat4  projection = glm::ortho( -aspect, aspect, -1.0f, 1.0f, -1.0f, 1.0f );
	// glUniformMatrix4fv( Projection, 1, GL_FALSE, glm::value_ptr(projection) );
	vp_width = width;
	vp_height = height;
	glUniform2f( Window, width, height );
	drawing_y = 0;
}
