#include <Windows.h>
#include <iostream>
#include <GL/glew.h>
#include <GL/GL.h>
#include <GL/freeglut.h>

#define GLFW_INCLUDE_GLU
#define GLFW_DLL
#include <GLFW/glfw3.h>
#include <vector>

#define GLM_SWIZZLE
#include <glm/glm.hpp>
#include <glm/gtc/constants.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtx/string_cast.hpp>

using namespace glm;

// -------------------------------------------------
// Global Variables
// -------------------------------------------------
int Width = 512;
int Height = 512;
std::vector<float> OutputImage;
// -------------------------------------------------

class Ray
{
public:
	vec3 p;
	vec3 d;
	Ray(vec3 p, vec3 d) {
		this->p = p;
		this->d = normalize(d);
	}

	vec3 getPos(float t) {
		return p + t * d;
	}
};

class Light
{
public:
	vec3 p;
	Light(vec3 p) {
		this->p = p;
	}
};

class Camera
{
public:
	vec3 e;
	vec3 u, v, w;
	float l, r, b, t, d;

	Camera(vec3 e) {
		this->e = e;
		u = vec3(1, 0, 0); v = vec3(0, 1, 0); w = vec3(0, 0, 1);
		l = -0.1f; r = 0.1f; b = -0.1f; t = 0.1f; d = 0.1f;
	}

	Ray* getRay(float ix, float iy) const {
		float uPos = l + (r - l) * (ix + 0.5f) / Width;
		float vPos = b + (t - b) * (iy + 0.5f) / Height;
		vec3 dir = normalize(-d * w + uPos * u + vPos * v);
		return new Ray(e, dir);
	}
};

struct RayHit
{
public:
	void* hitSurf;
	float t;

	RayHit(void* hitSurf, float t) {
		this->hitSurf = hitSurf;
		this->t = t;
	}
};

class Surface
{
protected:
	vec3 ka;
	vec3 kd;
	vec3 ks;
	float sp;
public:
	virtual vec3 getNormal(vec3 p) = 0;
	virtual RayHit* intersect(Ray* ray, float tMin, float tMax) = 0;
	virtual vec3 shade(Ray* ray, vec3 point, vec3 normal, Light* light, std::vector<Surface*>* surfs) = 0;
};

bool traceShadowRay(const vec3& start, const vec3& lightDir, float maxDist, std::vector<Surface*>* surfs) {
	Ray* shadowRay = new Ray(start, normalize(lightDir));
	bool inShadow = false;

	for (Surface* surf : *surfs) {
		float t;
		RayHit* rayHit = surf->intersect(shadowRay, 0.0f, maxDist);
		if (rayHit->hitSurf != nullptr) {
			inShadow = true;
			delete rayHit;
			break;
		}
		delete rayHit;
	}
	delete shadowRay;
	return inShadow;
}

vec3 gammaCorrection(vec3 color, float gamma) {
	return pow(color, vec3(1.0 / gamma));
}

class Plane : public Surface
{
	vec3 n;
	vec3 p;
public:
	Plane(vec3 n, vec3 p, vec3 ka, vec3 kd, vec3 ks, float sp) {
		this->n = normalize(n);
		this->p = p;
		this->ka = ka;
		this->kd = kd;
		this->ks = ks;
		this->sp = sp;
	}

	vec3 getNormal(vec3 p) override {
		return normalize(n);
	}

	RayHit* intersect(Ray* ray, float tMin, float tMax) override {
		float denom = dot(n, ray->d);
		float t = dot(p - ray->p, n) / denom;
		if (t >= tMin && t <= tMax) {
			RayHit* rayHit = new RayHit(this, t);
			return rayHit;
		}
		return new RayHit(nullptr, 0);
	}

	vec3 shade(Ray* ray, vec3 point, vec3 normal, Light* light, std::vector<Surface*>* surfs) override {
		// Ambient
		vec3 result = this->ka;

		vec3 l = normalize(light->p - point);
		bool inShadow = traceShadowRay(point + normal * 0.001f, l, length(light->p - point), surfs);

		if (!inShadow) {
			vec3 v = normalize(-ray->d);
			vec3 l = normalize(light->p - point);

			// Diffuse
			vec3 diffuse = this->kd * max(dot(normal, l), 0.0f);

			// Specular
			vec3 hd = normalize(v + l);
			float spec = pow(max(dot(normal, hd), 0.0f), this->sp);
			vec3 specular = this->ks * spec;

			result += diffuse + specular;
		}
		return gammaCorrection(result, 2.2);
	}
};

class Sphere : public Surface
{
	vec3 center;
	float radius;

public:
	Sphere(vec3 center, float radius, vec3 ka, vec3 kd, vec3 ks, float sp) {
		this->center = center;
		this->radius = radius;
		this->ka = ka;
		this->kd = kd;
		this->ks = ks;
		this->sp = sp;
	}

	vec3 getNormal(vec3 p) override {
		return normalize(p - center);
	}

	RayHit* intersect(Ray* ray, float tMin, float tMax) override {
		vec3 oc = ray->p - center;
		float b = -dot(ray->d, oc);
		float c = dot(oc, oc) - radius * radius;
		float D = b * b - c;
		if (D >= 0) {
			float rD = sqrtf(D);
			float t = min(b + rD, b - rD);
			if ((t >= tMin && t <= tMax)) {
				Surface* hitSurf = this;
				RayHit* rayHit = new RayHit(hitSurf, t);
				return rayHit;
			}
		}

		RayHit* rayHit = new RayHit(nullptr, 0);
 		return rayHit;
	}

	vec3 shade(Ray* ray, vec3 point, vec3 normal, Light* light, std::vector<Surface*>* surfs) override {
		// Ambient
		vec3 result = this->ka;

		vec3 l = normalize(light->p - point);
		bool inShadow = traceShadowRay(point + normal * 0.001f, l, length(light->p - point), surfs);

		if (!inShadow) {
			vec3 v = normalize(-ray->d);
			vec3 l = normalize(light->p - point);

			// Diffuse
			vec3 diffuse = this->kd * max(dot(normal, l), 0.0f);

			// Specular
			vec3 hd = normalize(v + l);
			float spec = pow(max(dot(normal, hd), 0.0f), this->sp);
			vec3 specular = this->ks * spec;

			result += diffuse + specular;
		}
		return gammaCorrection(result, 2.2);
	}
};

class Scene {
public:
	std::vector<Surface*> surfs;
	Light* light = new Light(vec3(-4.0f, 4.0f, -3.0f));

	vec3 trace(Ray* ray, float tMin, float tMax) {
		vec3 c = vec3(0, 0, 0);
		for (Surface* s : surfs) {
			float t;
			RayHit* rayHit = s->intersect(ray, tMin, tMax);
			if (rayHit->hitSurf != nullptr) {
				vec3 hitPoint = ray->getPos(rayHit->t);
				vec3 normal = s->getNormal(hitPoint);
				c = s->shade(ray, hitPoint, normal, light, &surfs);
			}
			delete rayHit;
		}
		return c;
	}
};

const int NUM_SAMPLES = 64;

vec3 getColorWithAntiAliasing(int x, int y, Camera* camera, Scene* scene) {
	vec3 color(0.0f);

	for (int i = 0; i < NUM_SAMPLES; i++) {
		float dx = float(rand()) / RAND_MAX;
		float dy = float(rand()) / RAND_MAX;

		float sampleX = x + dx;
		float sampleY = y + dy;

		Ray* ray = camera->getRay(sampleX, sampleY);
		vec3 sampleColor = scene->trace(ray, 0.0f, INFINITY);
		delete ray;

		color += sampleColor;
	}
	return color / float(NUM_SAMPLES);
}


void render()
{
	//Create our image. We don't want to do this in 
	//the main loop since this may be too slow and we 
	//want a responsive display of our beautiful image.
	//Instead we draw to another buffer and copy this to the 
	//framebuffer using glDrawPixels(...) every refresh
	OutputImage.clear();

	Scene scene;

	Plane p = Plane(vec3(0, 1, 0), vec3(0, -2, 0), vec3(0.2f, 0.2f, 0.2f), vec3(1.0f, 1.0f, 1.0f), vec3(0.0f, 0.0f, 0.0f), 0.0f);
	Sphere s1 = Sphere(vec3(-4.0f, 0.0f, -7.0f), 1.0f, vec3(0.2f, 0.0f, 0.0f), vec3(1.0f, 0.0f, 0.0f), vec3(0.0f, 0.0f, 0.0f), 0.0f);
	Sphere s2 = Sphere(vec3(0.0f, 0.0f, -7.0f), 2.0f, vec3(0.0f, 0.2f, 0.0f), vec3(0.0f, 0.5f, 0.0f), vec3(0.5f, 0.5f, 0.5f), 32.0f);
	Sphere s3 = Sphere(vec3(4.0f, 0.0f, -7.0f), 1.0f, vec3(0.0f, 0.0f, 0.2f), vec3(0.0f, 0.0f, 1.0f), vec3(0.0f, 0.0f, 0.0f), 0.0f);
	scene.surfs.push_back(&p);
	scene.surfs.push_back(&s1);
	scene.surfs.push_back(&s2);
	scene.surfs.push_back(&s3);

	Camera c = Camera(vec3(0, 0, 0));

	for (int j = 0; j < Height; ++j) 
	{
		for (int i = 0; i < Width; ++i) 
		{
			// ---------------------------------------------------
			// --- Implement your code here to generate the image
			// ---------------------------------------------------

			vec3 color = glm::vec3(0.5f, 0.5f, 0.5f); // grey color [0,1] in RGB channel

			//color = scene.trace(c.getRay(i, j), 0.0f, INFINITY);

			color = getColorWithAntiAliasing(i, j, &c, &scene);
			
			// set the color
			OutputImage.push_back(color.x); // R
			OutputImage.push_back(color.y); // G
			OutputImage.push_back(color.z); // B
		}
	}
}


void resize_callback(GLFWwindow*, int nw, int nh) 
{
	//This is called in response to the window resizing.
	//The new width and height are passed in so we make 
	//any necessary changes:
	Width = nw;
	Height = nh;
	//Tell the viewport to use all of our screen estate
	glViewport(0, 0, nw, nh);

	//This is not necessary, we're just working in 2d so
	//why not let our spaces reflect it?
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();

	glOrtho(0.0, static_cast<double>(Width)
		, 0.0, static_cast<double>(Height)
		, 1.0, -1.0);

	//Reserve memory for our render so that we don't do 
	//excessive allocations and render the image
	OutputImage.reserve(Width * Height * 3);
	render();
}


int main(int argc, char* argv[])
{
	// -------------------------------------------------
	// Initialize Window
	// -------------------------------------------------

	GLFWwindow* window;

	/* Initialize the library */
	if (!glfwInit())
		return -1;

	/* Create a windowed mode window and its OpenGL context */
	window = glfwCreateWindow(Width, Height, "OpenGL Viewer", NULL, NULL);
	if (!window)
	{
		glfwTerminate();
		return -1;
	}

	/* Make the window's context current */
	glfwMakeContextCurrent(window);

	//We have an opengl context now. Everything from here on out 
	//is just managing our window or opengl directly.

	//Tell the opengl state machine we don't want it to make 
	//any assumptions about how pixels are aligned in memory 
	//during transfers between host and device (like glDrawPixels(...) )
	glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
	glPixelStorei(GL_PACK_ALIGNMENT, 1);

	//We call our resize function once to set everything up initially
	//after registering it as a callback with glfw
	glfwSetFramebufferSizeCallback(window, resize_callback);
	resize_callback(NULL, Width, Height);

	/* Loop until the user closes the window */
	while (!glfwWindowShouldClose(window))
	{
		//Clear the screen
		glClear(GL_COLOR_BUFFER_BIT);

		// -------------------------------------------------------------
		//Rendering begins!
		glDrawPixels(Width, Height, GL_RGB, GL_FLOAT, &OutputImage[0]);
		//and ends.
		// -------------------------------------------------------------

		/* Swap front and back buffers */
		glfwSwapBuffers(window);

		/* Poll for and process events */
		glfwPollEvents();

		//Close when the user hits 'q' or escape
		if (glfwGetKey(window, GLFW_KEY_ESCAPE) == GLFW_PRESS
			|| glfwGetKey(window, GLFW_KEY_Q) == GLFW_PRESS)
		{
			glfwSetWindowShouldClose(window, GL_TRUE);
		}
	}

	glfwDestroyWindow(window);
	glfwTerminate();
	return 0;
}
