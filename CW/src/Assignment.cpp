#include <CanvasTriangle.h>
#include <DrawingWindow.h>
#include <Utils.h>
#include <fstream>
#include <vector>
#include <glm/glm.hpp>
#include <CanvasPoint.h>
#include <CanvasTriangle.h>
#include <Colour.h>
#include <ModelTriangle.h>
#include <TextureMap.h>
#include <TexturePoint.h>
#include <Utils.h>
#include <unordered_map>
#include <RayTriangleIntersection.h>
#include <math.h>


#define WIDTH 640
#define HEIGHT 480

bool orbiting = false;

glm::vec3 camPos(0.0, 0.0, 4.0);
glm::vec3 light(1.0, 1.0, 2.0);


std::vector<glm::vec3> lights{};


// for sphere
// glm::vec3 camPos(0.0, 1.3, 3.5);
// glm::vec3 light(1.0, 2.0, 2.8);

glm::mat3 camOrientation(
	glm::vec3(1.0,0.0,0.0),
	glm::vec3(0.0,1.0,0.0),
	glm::vec3(0.0,0.0,1.0)
);


int planeMultiplier = 200;

//lighting settings
bool proximity = true;
bool incidence = true;
bool specular = true;


void initialiseLights(int size){
	glm::vec3 light(0.0,1.0,0.5);
	lights.push_back(light);
	for (int i=0; i<size; i++) {
		float j = (i + 1) * 0.025;
		lights.push_back(light + glm::vec3(j, j, j));
		lights.push_back(light + glm::vec3(j, j, -j));
		lights.push_back(light + glm::vec3(j, -j, j));
		lights.push_back(light + glm::vec3(-j, j, j));
		lights.push_back(light + glm::vec3(j, -j, -j));
		lights.push_back(light + glm::vec3(-j, j, -j));
		lights.push_back(light + glm::vec3(-j, -j, j));
		lights.push_back(light + glm::vec3(-j, -j, -j));
	}
}

std::unordered_map<std::string, Colour> loadMtlFile(const std::string &filename, std::unordered_map<std::string, TextureMap> &textures) {
	std::unordered_map<std::string, Colour> colours;

	std::ifstream inputStream(filename, std::ifstream::in);
	std::string nextLine;
	std::string colour_name;

	while (std::getline(inputStream, nextLine)) {
		std::vector<std::string> line = split(nextLine, ' ');

		if (line[0] == "newmtl") {
			colour_name = line[1];	
		} else if (line[0] == "Kd") {
			Colour colour(int(std::stof(line[1])*255),int(std::stof(line[2])*255),int(std::stof(line[3])*255));
			colours.insert({colour_name, colour});
		} else if(line[0] == "map_Kd") {
			Colour colour = colours[colour_name];
			colour.name = line[1];
			textures.insert({line[1], TextureMap(line[1])});
			colours[colour_name] = colour;
		}
	}
	inputStream.close();

	return colours;
}

std::vector<ModelTriangle> vertexNormals(std::vector<ModelTriangle> triangles) {
	for(int i = 0; i < triangles.size(); i++) {
		ModelTriangle triangle = triangles[i];
		std::vector<glm::vec3> vertex_normals;
		for(int v = 0; v < triangle.vertices.size(); v++) {
			glm::vec3 vertex = triangle.normal;
			int count = 1;
			for(int j = 0; j < triangles.size(); j++) {
				ModelTriangle tri = triangles[j];
				for(int u = 0; u < tri.vertices.size(); u++) {
					if(i != j && triangle.vertices[v].x == tri.vertices[u].x && triangle.vertices[v].y == tri.vertices[u].y && triangle.vertices[v].z == tri.vertices[u].z) {
						if (std::acos(glm::dot(triangle.normal, tri.normal) / (length(triangle.normal) * length(tri.normal))) < M_PI / 4) { //??
							vertex += tri.normal;
							count ++;
						}
					}
				}
			}
			vertex = vertex / float(count);
			triangles[i].vertexNormals[v] = normalize(vertex);
		}
	}
	return triangles;
}


std::vector<ModelTriangle> loadObjFile(const std::string &filename, float scale, std::unordered_map<std::string, TextureMap> &textures) {
	std::vector<glm::vec3> vertices;
	std::vector<ModelTriangle> faces;
	std::vector<TexturePoint> texturePoints;
	std::unordered_map<std::string, Colour> materials;
	std::vector<glm::vec3> normals;
	bool mirror = false;

	std::ifstream inputStr(filename, std::ifstream::in);
	std::string nextLine;
	Colour colour;
	while (std::getline(inputStr, nextLine)) { //extracts from inputStr and stores into nextLine
		std::vector<std::string> vector = split(nextLine, ' '); //split line by spaces
		if (vector[0] == "mtllib") {
			materials = loadMtlFile(vector[1], textures);
		}else if (vector[0] == "usemtl") {
			colour = materials[vector[1]];
			if(vector[1] == "Mirror") mirror = true;
		}else if (vector[0] == "v") {
			vertices.push_back(glm::vec3(
				std::stof(vector[1]) * scale, //string to float
				std::stof(vector[2]) * scale,
				std::stof(vector[3]) * scale
			));
		}else if (vector[0] == "vn") {
			normals.push_back(glm::vec3(
				std::stof(vector[1]),
				std::stof(vector[2]),
				std::stof(vector[3])
			));
		}else if (vector[0] == "f") { //indexed from 1
			ModelTriangle triangle = ModelTriangle();
			for (int i=0; i < 3; i++) {
				std::vector<std::string> v = split(vector[i + 1], '/');
				triangle.vertices[i] = vertices[std::stoi(v[0]) - 1];
				if (v.size() > 1 && v[1] != "") triangle.texturePoints[i] = texturePoints[std::stoi(v[1]) - 1];
				if (v.size() > 2 && v[2] != "") triangle.vertexNormals[i] = normals[std::stoi(v[2]) - 1];

			}
			triangle.colour = colour;
			triangle.normal = glm::normalize(glm::cross(glm::vec3(triangle.vertices[1] - triangle.vertices[0]), glm::vec3(triangle.vertices[2] - triangle.vertices[0])));
			triangle.mirror = mirror;
			faces.push_back(triangle);
		}else if(vector[0] == "vt") {
			texturePoints.push_back(TexturePoint(stof(vector[1]), stof(vector[2])));
		}
	}

	if(normals.empty()) {
		faces = vertexNormals(faces);
	} 
	inputStr.close();
	return faces;
}

std :: vector<float> interpolateSingleFloats(float from, float to, int numberOfValues) {
	float step = (to - from) / (numberOfValues -1);
	std :: vector<float> v;
	for (int i = 0; i < numberOfValues; i++){
		v.push_back(from + i*step);
	}
	return v;
}

std::vector<CanvasPoint> interpolateRoundPoints(CanvasPoint from, CanvasPoint to, int numberOfValues) {
	std::vector<CanvasPoint> interpolatedValues;
	std::vector<float> xs = interpolateSingleFloats(from.x, to.x, numberOfValues);
	std::vector<float> ys = interpolateSingleFloats(from.y, to.y, numberOfValues);
	std::vector<float> depths = interpolateSingleFloats(from.depth, to.depth, numberOfValues);
	for (int i=0; i<numberOfValues; i++) {
		float x = i==0 ? floor(xs[i]) : i==numberOfValues-1 ? ceil(xs[i]) : round(xs[i]);
		float y = i==0 ? floor(ys[i]) : i==numberOfValues-1 ? ceil(ys[i]) : round(ys[i]);
		interpolatedValues.push_back(CanvasPoint(x, y, depths[i]));
	}
	return interpolatedValues;
}

std::vector<TexturePoint> interpolateRoundPoints(TexturePoint from, TexturePoint to, int numberOfValues) {
	std::vector<TexturePoint> interpolatedValues;
	std::vector<float> xs = interpolateSingleFloats(from.x, to.x, numberOfValues);
	std::vector<float> ys = interpolateSingleFloats(from.y, to.y, numberOfValues);
	
	for (int i=0; i<numberOfValues; i++) {
		float x = i==0 ? floor(xs[i]) : i==numberOfValues-1 ? ceil(xs[i]) : round(xs[i]);
		float y =i==0 ? floor(ys[i]) : i==numberOfValues-1 ? ceil(ys[i]) : round(ys[i]);
		interpolatedValues.push_back(TexturePoint(x, y));
	}
	return interpolatedValues;
}

void sortTriangle(CanvasTriangle &triangle){
	if (triangle[0].y > triangle[1].y) {
		std::swap(triangle[0], triangle[1]);
	}
	if (triangle[1].y > triangle[2].y) {
		std::swap(triangle[1], triangle[2]);
		if (triangle[0].y > triangle[1].y) {
			std::swap(triangle[0], triangle[1]);
		}	
	}
}

void drawLine(DrawingWindow &window, CanvasPoint from, CanvasPoint to, Colour colour, std::vector<std::vector<float>> &depthBuffer) {
	float numberOfSteps = fmax(fmax(abs(to.x - from.x), abs(to.y - from.y)), 1.0f);
	std :: vector<CanvasPoint> points = interpolateRoundPoints(from, to, numberOfSteps + 1);
	for (int i=0; i<numberOfSteps; i++) {
		int x = points[i].x;
		int y = points[i].y;
		if(x >= 0 && x < window.width && y >= 0 && y < window.height) {
			float pointDepth = 1 / -points[i].depth;
			if (pointDepth > depthBuffer[y][x]) {
				depthBuffer[y][x] = pointDepth;
				window.setPixelColour(x, y, (255 << 24) + (colour.red << 16) + (colour.green << 8) + colour.blue);
			}
		}
	}
}


void drawStrokedTriangle(DrawingWindow &window, CanvasTriangle triangle, Colour colour, std::vector<std::vector<float>> &depthBuffer) {
	drawLine(window, triangle[0], triangle[1], colour, depthBuffer);
	drawLine(window, triangle[0], triangle[2], colour, depthBuffer);
	drawLine(window, triangle[1], triangle[2], colour, depthBuffer);
}

void drawFilledTriangle(DrawingWindow &window, CanvasTriangle triangle, Colour colour, std::vector<std::vector<float>> &depthBuffer) {
	sortTriangle(triangle);
	std :: vector<CanvasPoint> start = interpolateRoundPoints(triangle[0], triangle[1], triangle[1].y - triangle[0].y + 1);
	if (triangle[2].y - triangle[1].y + 1 > 1){
		start.pop_back(); // Last row duplicated by next triangle so pop
		std :: vector<CanvasPoint> start2 = interpolateRoundPoints(triangle[1], triangle[2], triangle[2].y - triangle[1].y + 1);
		start.insert(start.end(), start2.begin(), start2.end());
	}
	
	std :: vector<CanvasPoint> end = interpolateRoundPoints(triangle[0], triangle[2], triangle[2].y - triangle[0].y + 1);

	// Draw the filled in triangle
	for (int i=0; i<=triangle[2].y - triangle[0].y; i++) {
		drawLine(window, start[i], end[i], colour, depthBuffer);
	}

	drawStrokedTriangle(window, triangle, colour, depthBuffer); //outline

}


void drawTexturedTriangle(DrawingWindow &window, CanvasTriangle triangle, TextureMap texMap, std::vector<std::vector<float>> &depthBuffer) {

	sortTriangle(triangle);
	std :: vector<CanvasPoint> xStart = interpolateRoundPoints(triangle[0], triangle[1], triangle[1].y - triangle[0].y + 1);
	if (triangle[2].y - triangle[1].y + 1 > 1){
		xStart.pop_back();
		std :: vector<CanvasPoint> xStart2 = interpolateRoundPoints(triangle[1], triangle[2], triangle[2].y - triangle[1].y + 1);
		xStart.insert(xStart.end(), xStart2.begin(), xStart2.end());
	}
	
	std :: vector<CanvasPoint> xEnd = interpolateRoundPoints(triangle[0], triangle[2], triangle[2].y - triangle[0].y + 1);

	std :: vector<TexturePoint> textureStarts = interpolateRoundPoints(triangle[0].texturePoint, triangle[1].texturePoint, triangle[1].y - triangle[0].y + 1);
	if (triangle[2].y - triangle[1].y + 1 > 1){
		textureStarts.pop_back();
		std :: vector<TexturePoint> textureStarts2 = interpolateRoundPoints(triangle[1].texturePoint, triangle[2].texturePoint, triangle[2].y - triangle[1].y + 1);
		textureStarts.insert(textureStarts.end(), textureStarts2.begin(), textureStarts2.end());
	}
	std :: vector<TexturePoint> textureEnds = interpolateRoundPoints(triangle[0].texturePoint, triangle[2].texturePoint, triangle[2].y - triangle[0].y + 1);

	for(int i=0; i <= triangle[2].y - triangle[0].y; i++) {
		float numberOfSteps = fmax(abs(xStart[i].x - xEnd[i].x), 1.0f);
		//interpolate between start and end of rake to find texture point for pixel
		std :: vector<TexturePoint> texPoints = interpolateRoundPoints(textureStarts[i], textureEnds[i], numberOfSteps + 1);
		std :: vector<CanvasPoint> points = interpolateRoundPoints(xStart[i], xEnd[i], numberOfSteps + 1);

		for (float j=0.0; j <= numberOfSteps; j++) {
			int x = round(points[j].x);
			int y = round(points[j].y);
			if(x >= 0 && x < window.width && y >= 0 && y < window.height) {
			
				float pointDepth = 1 / -points[j].depth;

				if (pointDepth > depthBuffer[y][x]) {
					depthBuffer[y][x] = pointDepth;
					uint32_t col = texMap.pixels[round(texPoints[j].y) * texMap.width + round(texPoints[j].x)];
					window.setPixelColour(x, y, col);
				}
			}
		}
	}

}



CanvasPoint getCanvasIntersectionPoint(DrawingWindow &window, glm::vec3 vertexPosition, float focalLength){
	glm::vec3 vertex = (vertexPosition - camPos) * camOrientation; //camera co-ordinate system
	float u = -round(planeMultiplier*focalLength * (vertex.x / vertex.z)) + (window.width / 2);
	float v = round(planeMultiplier*focalLength * (vertex.y / vertex.z)) + (window.height / 2);
	float z = vertex.z;

	return CanvasPoint(u,v,z);
}

glm::mat3 rotateY(float angle) {
	return glm::mat3(
		cos(angle), 0.0, -sin(angle),
		0.0, 1.0, 0.0,
		sin(angle), 0.0, cos(angle)
	);
}

glm::mat3 rotateX(float angle) {
	return glm::mat3(
		1.0, 0.0, 0.0,
		0.0, cos(angle), sin(angle),
		0.0, -sin(angle), cos(angle)
	);
}	



void lookAt() {
	glm::vec3 lookAtPoint(0, 0, 0);
	glm::vec3 forward = glm::normalize(camPos - lookAtPoint);
	glm::vec3 right = glm::normalize(glm::cross(glm::vec3(0, 1, 0), forward));
	glm::vec3 up = glm::normalize(glm::cross(forward, right));
	camOrientation =  glm::mat3(right, up, forward);
}


void orbit(bool orb){ 
	if (orb){
		camPos =  rotateY(-0.1) * camPos;
		lookAt();
	}	
}

void reset_camera() {
	camPos = glm::vec3(0.0,0.0,4.0);
	camOrientation = glm::mat3(glm::vec3(1.0,0.0,0.0),glm::vec3(0.0,1.0,0.0),glm::vec3(0.0,0.0,1.0));
}

bool is_shadow(RayTriangleIntersection intersect, std::vector<ModelTriangle> triangles, glm::vec3 light) {
	glm::vec3 shadowRay = light-intersect.intersectionPoint;
	for(int i = 0; i < triangles.size(); i++) {
		if (i != intersect.triangleIndex){
			ModelTriangle triangle = triangles[i];
			glm::vec3 e0 = triangle.vertices[1] - triangle.vertices[0];
			glm::vec3 e1 = triangle.vertices[2] - triangle.vertices[0];
			glm::vec3 SPVector = intersect.intersectionPoint - triangle.vertices[0];
			glm::mat3 DEMatrix(-glm::normalize(shadowRay), e0, e1);
			glm::vec3 possibleSolution = glm::inverse(DEMatrix) * SPVector;
			float t = possibleSolution.x, u = possibleSolution.y, v = possibleSolution.z;

			if((u >= 0.0) && (u <= 1.0) && (v >= 0.0) && (v <= 1.0) && (u + v) <= 1.0) {
				if(t < glm::length(shadowRay) && t > 0.05) {
					return true;
				}
			}
		}
	}
	return false;
}

RayTriangleIntersection getClosestReflection(glm::vec3 inter, glm::vec3 direction, std::vector<ModelTriangle> triangles, int index) {
	RayTriangleIntersection intersection;
	intersection.distanceFromCamera = std::numeric_limits<float>::infinity();

	for(int i = 0; i < triangles.size(); i++) {
		if (i != index){
			ModelTriangle triangle = triangles[i];
			glm::vec3 e0 = triangle.vertices[1] - triangle.vertices[0];
			glm::vec3 e1 = triangle.vertices[2] - triangle.vertices[0];
			glm::vec3 SPVector = inter - triangle.vertices[0];
			glm::mat3 DEMatrix(-direction, e0, e1);
			glm::vec3 possibleSolution = glm::inverse(DEMatrix) * SPVector;
			float t = possibleSolution.x, u = possibleSolution.y, v = possibleSolution.z;

			if((u >= 0.0) && (u <= 1.0) && (v >= 0.0) && (v <= 1.0) && (u + v) <= 1.0 && intersection.distanceFromCamera > t && t > 0.0) {
				glm::vec3 intersect = triangle.vertices[0]+u*e0+v*e1;
				intersection =  RayTriangleIntersection(intersect, t, triangle, i);
				intersection.u = u;
				intersection.v = v;
			}
		}
	}
	return intersection;
}


RayTriangleIntersection getClosestIntersection(glm::vec3 rayDirection, std::vector<ModelTriangle> triangles) {
	RayTriangleIntersection intersection;
	intersection.distanceFromCamera = std::numeric_limits<float>::infinity();
	float distance = std::numeric_limits<float>::infinity();

	for(int i = 0; i < triangles.size(); i++) {
		ModelTriangle triangle = triangles[i];
		glm::vec3 e0 = triangle.vertices[1] - triangle.vertices[0];
		glm::vec3 e1 = triangle.vertices[2] - triangle.vertices[0];
		glm::vec3 SPVector = camPos - triangle.vertices[0];
		glm::mat3 DEMatrix(-rayDirection, e0, e1);
		glm::vec3 possibleSolution = glm::inverse(DEMatrix) * SPVector;
		float t = possibleSolution.x, u = possibleSolution.y, v = possibleSolution.z;
		if((u >=  0.0) && (u <= 1.0) && (v >= 0.0) && (v <= 1.0) && (u + v) <= 1.0 && t < distance && t > 0.0) {
			glm::vec3 intersect = triangle.vertices[0]+u*e0+v*e1;
			distance = t;
			if(triangle.mirror) {
				glm::vec3 normal = normalize(triangle.normal);
				glm::vec3 reflection_ray = normalize(rayDirection - (normal * 2.0f * glm::dot(rayDirection, normal)));
				intersection = getClosestReflection(intersect, reflection_ray, triangles, i);
				if(isinf(intersection.distanceFromCamera)) intersection.inf = true;
			} else {
				intersection =  RayTriangleIntersection(intersect, t, triangle, i);
				intersection.u = u;
				intersection.v = v;
			}
		}
	}
	return intersection;
}


float getBrightness(glm::vec3 intersectionPoint, glm::vec3 normal, glm::vec3 light) {
	float intensity = 30;
	float specularScale = 256;

	glm::vec3 lightRay = light - intersectionPoint;
	float length = glm::length(lightRay);

	glm::vec3 cameraRay = glm::normalize(camOrientation * (camPos - intersectionPoint));

	float brightness = proximity ? intensity/(4 * M_PI * length*length) : 1.0; // Proximity lighting
	float angleOfIncidence = incidence ? glm::dot(glm::normalize(lightRay), normal) : 1.0; // Angle of Incidence Lighting
	brightness *= std::max(0.0f, angleOfIncidence);

	glm::vec3 angleOfReflection = glm::normalize(lightRay) - ((2.0f*normal)*glm::dot(glm::normalize(lightRay), normal));
	float specularLight = specular ? std::pow(glm::dot(glm::normalize(angleOfReflection), cameraRay), specularScale) : 0.0; // Specular Lighting

	brightness += std::max(0.0f, specularLight*0.2f);
	brightness = std::max(0.2f, std::min(1.0f, brightness)); // Ambient Lighting

	return brightness;
}

float gouraud(RayTriangleIntersection intersection, glm::vec3 light) {
	float intensity = 30;
	float specularScale = 256;

	glm::vec3 lightRay = light - intersection.intersectionPoint;
	float length = glm::length(lightRay);
	glm::vec3 cameraRay = glm::normalize(camOrientation * (camPos - intersection.intersectionPoint));

	ModelTriangle triangle = intersection.intersectedTriangle;
	
	std::vector<glm::vec3> reflections;
	std::vector<float> incidences;

	float brightness = proximity ? intensity/(4 * M_PI * length*length) : 1.0; //proximity lighting

	for(int i = 0; i < triangle.vertexNormals.size(); i++) {
		incidences.push_back(glm::dot(triangle.vertexNormals[i], glm::normalize(lightRay)));
		reflections.push_back(glm::normalize(glm::normalize(lightRay) - ((2.0f*triangle.vertexNormals[i])*glm::dot(glm::normalize(lightRay), triangle.vertexNormals[i]))));
	}
	float angleOfIncidence = incidence ? (1 - intersection.u - intersection.v) * incidences[0] + intersection.u * incidences[1] + intersection.v * incidences[2] : 1.0;
	brightness *= std::max(0.0f, angleOfIncidence); 

	glm::vec3 angleOfReflection = (1 - intersection.u - intersection.v) * reflections[0] + intersection.u * reflections[1] + intersection.v * reflections[2];
	float specularLighting = specular ? std::pow(glm::dot(angleOfReflection, cameraRay), specularScale) : 0.0;

	brightness += std::max(0.0f, specularLighting);
	brightness = std::max(0.2f, std::min(1.0f, brightness));

	
	return brightness;
}

float phong(RayTriangleIntersection intersection, glm::vec3 light) {
	float specularScale = 256;
	float intensity = 30;

	glm::vec3 lightRay = light - intersection.intersectionPoint;
	float length = glm::length(lightRay);
	glm::vec3 cameraRay = glm::normalize(camOrientation * (camPos - intersection.intersectionPoint));


	ModelTriangle triangle = intersection.intersectedTriangle;
	//interpolated normals
	glm::vec3 normal = (1 - intersection.u - intersection.v) * triangle.vertexNormals[0] + intersection.u * triangle.vertexNormals[1] + intersection.v * triangle.vertexNormals[2];
	
	float brightness = proximity ? intensity/(4 * M_PI * length*length) : 1.0; // Proximity lighting
	float angleOfIncidence =  incidence ? glm::dot(glm::normalize(lightRay),glm::normalize(normal)) : 1.0;
	brightness *= std::max(0.0f, angleOfIncidence); 


	glm::vec3 angleOfReflection = glm::normalize(glm::normalize(lightRay) - (glm::normalize(normal)*2.0f*glm::dot(glm::normalize(lightRay), glm::normalize(normal))));
	float specularLighting = specular ? std::pow(glm::dot(angleOfReflection, cameraRay), specularScale) : 0.0; // Specular Lighting

	brightness += std::max(0.0f, specularLighting);
	brightness = std::max(0.2f, std::min(1.0f, brightness)); //ambient
	
	return brightness;
}

void drawRayTrace(DrawingWindow &window, std::vector<ModelTriangle> triangles, float focalLength, int lightMode, std::unordered_map<std::string, TextureMap> textures) {
	window.clearPixels();
	orbit(orbiting); 
	for(int x = 0; x < window.width; x++) {
		for(int y = 0; y < window.height; y++) {
			glm::vec3 direction(x-float(window.width / 2), float(window.height / 2) - y, -focalLength*planeMultiplier);;
			glm::vec3 ray = normalize(camOrientation * (direction-camPos));
			RayTriangleIntersection intersection = getClosestIntersection(ray, triangles);

			Colour colour;

			if(!isinf(intersection.distanceFromCamera)) {
				float brightness = 0;
				for (int i = 0; i < lights.size(); i++) {
					if(!is_shadow(intersection, triangles, lights[i])){
						if (lightMode == 0){
							glm::vec3 normal = intersection.intersectedTriangle.normal;
							brightness += getBrightness(intersection.intersectionPoint, normal, lights[i]);
						}else if (lightMode == 1){
							brightness += gouraud(intersection, lights[i]);
						}else {
							brightness += phong(intersection, lights[i]);
						}
					}else{
						brightness += 0.15; //shadow brightness
					}
				}
				brightness /= lights.size();
				
				if (intersection.intersectedTriangle.colour.name != "") { //texture
					ModelTriangle triangle = intersection.intersectedTriangle; //triangles[intersection.triangleIndex];
                    TextureMap texture = textures[triangle.colour.name];
                    // interpolate the ratios of texture points
                    float ratioX = (1 - intersection.u - intersection.v) * triangle.texturePoints[0].x + intersection.u * triangle.texturePoints[1].x + intersection.v * triangle.texturePoints[2].x;
                    float ratioY = (1 - intersection.u - intersection.v) * triangle.texturePoints[0].y + intersection.u * triangle.texturePoints[1].y + intersection.v * triangle.texturePoints[2].y;
					float x = fmod(ratioX, 1) * texture.width; //get co-ordinates
					float y = texture.height - fmod(ratioY,1)*texture.height;
					uint32_t c = texture.pixels[round(y)*texture.width + round(x)];
					
                    uint8_t red = (c >> 16) & 0xff;
                    uint8_t green = (c >> 8) & 0xff;
                    uint8_t blue = c & 0xff;
					colour = Colour(red, green, blue);
				}else{
					colour = intersection.intersectedTriangle.colour;
				}

				if(intersection.inf) {
					colour = Colour(0, 0, 0);
				}

				colour.red *= brightness;
				colour.green *= brightness;
				colour.blue *= brightness;
				uint32_t c = (255 << 24) + (colour.red << 16) + (colour.green << 8) + colour.blue;
				window.setPixelColour(x, y, c);
				

			}
		}
	}
	
}


void drawRasterisedScene(DrawingWindow &window, std::vector<ModelTriangle> faces, float focalLength, int renderMode) {
	window.clearPixels();
	orbit(orbiting); 
	std::vector<std::vector<float>> depthBuffer (window.height, std::vector<float>(window.width, 0));

	for (int i=0; i<faces.size(); i++) {
		ModelTriangle face = faces[i];
		CanvasTriangle triangle = CanvasTriangle();
		for (int j=0; j<face.vertices.size(); j++) {
			glm::vec3 vertexPosition = face.vertices[j];
			triangle.vertices[j] = getCanvasIntersectionPoint(window,vertexPosition,focalLength);
		}
		
		if (renderMode == 0) {
			drawStrokedTriangle(window, triangle, Colour(255, 255, 255), depthBuffer); //wireframe
		}
		else if (face.colour.name != "") {
			

			TextureMap texture(face.colour.name);
			for(int j = 0; j < triangle.vertices.size(); j++) {
				triangle.vertices[j].texturePoint = face.texturePoints[j];
				triangle.vertices[j].texturePoint.x = fmod(triangle.vertices[j].texturePoint.x, 1 ) * texture.width;
				triangle.vertices[j].texturePoint.y = texture.height-fmod(triangle.vertices[j].texturePoint.y, 1) * texture.height;
			}
			drawTexturedTriangle(window, triangle, texture, depthBuffer);
		}
		else {
			drawFilledTriangle(window, triangle, faces[i].colour, depthBuffer);


		}
	}
}




void handleEvent(SDL_Event event, DrawingWindow &window, int &renderMode, int &lightMode) {
	if (event.type == SDL_KEYDOWN) {
		//camera movement
		if (event.key.keysym.sym == SDLK_e) camPos.y += 0.1; //up
		else if (event.key.keysym.sym == SDLK_q) camPos.y -= 0.1; //down
		else if (event.key.keysym.sym == SDLK_w) camPos.z -= 0.1; //forwards
		else if (event.key.keysym.sym == SDLK_s) camPos.z += 0.1; //backwards
		else if (event.key.keysym.sym == SDLK_d) camPos.x += 0.1; //right
		else if (event.key.keysym.sym == SDLK_a) camPos.x -= 0.1; //left
		// Order of transformation is matrix * vector
		//rotation
		else if (event.key.keysym.sym == SDLK_UP) {
			camPos = rotateY(0.1) * camPos; //rotate Y C
			lookAt();
		}else if (event.key.keysym.sym == SDLK_RIGHT){
			camPos = rotateX(0.1) * camPos; //rotate X C
			lookAt();
		}else if (event.key.keysym.sym == SDLK_DOWN){
			camPos = rotateY(-0.1) * camPos; //rotate Y antiC
			lookAt();
		}
		else if (event.key.keysym.sym == SDLK_LEFT){
			camPos = rotateX(-0.1) * camPos; //rotate X antiC
			lookAt();
		}//orientation
		else if (event.key.keysym.sym == SDLK_i){
			camOrientation = rotateY(0.1) * camOrientation ; //panning
		}else if (event.key.keysym.sym == SDLK_k){
			camOrientation = rotateY(-0.1) * camOrientation;
		}else if (event.key.keysym.sym == SDLK_l){
			camOrientation = rotateX(0.1) * camOrientation; //tilting
		}else if (event.key.keysym.sym == SDLK_j){
			camOrientation = rotateX(-0.1) * camOrientation;
		}

		//orbit
		else if (event.key.keysym.sym == SDLK_o) orbiting = (orbiting) ? false : true;
		else if (event.key.keysym.sym == SDLK_r) reset_camera();

		//render Mode
		else if (event.key.keysym.sym == SDLK_0) renderMode = 0;
		else if (event.key.keysym.sym == SDLK_1) renderMode = 1;
		else if (event.key.keysym.sym == SDLK_2) renderMode = 2;
		//light Mode
		else if (event.key.keysym.sym == SDLK_3) lightMode = 0; //normal
		else if (event.key.keysym.sym == SDLK_4) lightMode = 1; //gouraud
		else if (event.key.keysym.sym == SDLK_5) lightMode = 2; //phong
		else if (event.key.keysym.sym == SDLK_6) proximity = !proximity;
		else if (event.key.keysym.sym == SDLK_7) incidence = !incidence; 
		else if (event.key.keysym.sym == SDLK_8) specular = !specular;


	} else if (event.type == SDL_MOUSEBUTTONDOWN) {
		window.savePPM("output.ppm");
		window.saveBMP("output.bmp");
	}
}

int main(int argc, char *argv[]) {
	DrawingWindow window = DrawingWindow(WIDTH, HEIGHT, false);
	SDL_Event event;

	float vertexScale = 0.5;
	int renderMode = 1;
	int lightMode = 0;

	initialiseLights(3);

	std::unordered_map<std::string, TextureMap> textures;
	
	std::vector<ModelTriangle> triangles = loadObjFile("logo.obj", 0.003, textures);
	// std::vector<ModelTriangle> triangles = loadObjFile("textured-cornell-box.obj", vertexScale, textures);
	// std::vector<ModelTriangle> triangles = loadObjFile("cornell-box.obj", vertexScale, textures);



	 //std::vector<ModelTriangle> triangles = loadObjFile("sphere2.obj", vertexScale, textures);
	// triangles.insert(triangles.end(), sphere.begin(), sphere.end());
	// triangles.insert(triangles.end(), logo.begin(), logo.end());

	float focalLength = 2.0;
	while (true) {
		if (window.pollForInputEvents(event)){
			handleEvent(event, window, renderMode, lightMode);
		}
		if (renderMode == 2) drawRayTrace(window, triangles, focalLength, lightMode, textures);
		else drawRasterisedScene(window, triangles, focalLength, renderMode);

		window.renderFrame();
	}
}