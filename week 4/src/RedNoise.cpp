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

#define WIDTH 320
#define HEIGHT 240

std::unordered_map<std::string, Colour> loadMtlFile(const std::string &filename) {
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
		}
	}
	inputStream.close();
	return colours;
}

std::vector<ModelTriangle> loadObjFile(const std::string &filename, float scale, std::unordered_map<std::string, Colour> materials) {
	std::vector<glm::vec3> vertices;
	std::vector<ModelTriangle> faces;

	std::ifstream inputStr(filename, std::ifstream::in);
	std::string nextLine;
	Colour colour;
	while (std::getline(inputStr, nextLine)) { //extracts from inputStr and stores into nextLine
		std::vector<std::string> vector = split(nextLine, ' '); //split line by spaces
		if (vector[0] == "usemtl") {
			colour = materials[vector[1]];
		}else if (vector[0] == "v") {
			vertices.push_back(glm::vec3(
				std::stof(vector[1]) * scale, //string to float
				std::stof(vector[2]) * scale,
				std::stof(vector[3]) * scale
			));
		}
		else if (vector[0] == "f") { //indexed from 1
			faces.push_back(ModelTriangle(
				vertices[std::stoi(split(vector[1], '/')[0]) - 1],
				vertices[std::stoi(split(vector[2], '/')[0]) - 1],
				vertices[std::stoi(split(vector[3], '/')[0]) - 1],
				colour
			));
		}
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

template <typename T>
std::vector<CanvasPoint> interpolateRoundPoints(T from, T to, int numberOfValues) {
	std::vector<CanvasPoint> points;
	std :: vector<float> xs = interpolateSingleFloats(from.x, to.x, numberOfValues);
	std :: vector<float> ys = interpolateSingleFloats(from.y, to.y, numberOfValues);
	std :: vector<float> depths = interpolateSingleFloats(from.depth, to.depth, numberOfValues);
	for (int i=0; i<numberOfValues; i++) {
		points.push_back(CanvasPoint(round(xs[i]), round(ys[i]),  depths[i]));
	}
	return points;
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

void drawLine(DrawingWindow &window, CanvasPoint from, CanvasPoint to, Colour colour, std::vector<float> &depthBuffer) {
	float numberOfSteps = fmax(fmax(abs(to.x - from.x), abs(to.y - from.y)), 1);
	std :: vector<CanvasPoint> points = interpolateRoundPoints(from, to, numberOfSteps + 1);
	for (int i=0; i<=numberOfSteps; i++) {

		int depthIndex = (points[i].y * window.width) + points[i].x;
		float pointDepth = 1 / -points[i].depth;
		if (pointDepth > depthBuffer[depthIndex]) {
			depthBuffer[depthIndex] = pointDepth;
			window.setPixelColour(points[i].x, points[i].y, (255 << 24) + (colour.red << 16) + (colour.green << 8) + colour.blue);
		}
	}
}


void drawStrokedTriangle(DrawingWindow &window, CanvasTriangle triangle, Colour colour, std::vector<float> &depthBuffer) {
	drawLine(window, triangle[0], triangle[1], colour, depthBuffer);
	drawLine(window, triangle[0], triangle[2], colour, depthBuffer);
	drawLine(window, triangle[1], triangle[2], colour, depthBuffer);
}

void drawFilledTriangle(DrawingWindow &window, CanvasTriangle triangle, Colour colour, std::vector<float> &depthBuffer, bool outline) {
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

	if (outline) drawStrokedTriangle(window, triangle, Colour(255,255,255), depthBuffer); //outline

}


CanvasPoint getCanvasIntersectionPoint(DrawingWindow &window, glm::vec3 cameraPosition,glm::vec3 vertexPosition, float focalLength){

	int planeMultiplier = 450;
	glm::vec3 vertex = vertexPosition - cameraPosition;
	float u = -round(planeMultiplier*focalLength * (vertex.x / vertex.z)) + (window.width / 2);
	float v = round(planeMultiplier*focalLength * (vertex.y / vertex.z)) + (window.height / 2);
	float z = vertex.z;

	return CanvasPoint(u,v,z);
}

void pointcloud(DrawingWindow &window, glm::vec3 cameraPosition, std::vector<ModelTriangle> faces, float focalLength) {
	window.clearPixels();
	std::vector<float> depthBuffer = std::vector<float>(window.height * window.width, 0);

	for (int i=0; i<faces.size(); i++) {
		ModelTriangle face = faces[i];
		CanvasTriangle triangle = CanvasTriangle();
		for (int j=0; j<face.vertices.size(); j++) {
			glm::vec3 vertexPosition = face.vertices[j];
			triangle.vertices[j] = getCanvasIntersectionPoint(window,cameraPosition, vertexPosition, focalLength );
		}
		drawFilledTriangle(window, triangle, faces[i].colour, depthBuffer, false);
	}
	
}

void handleEvent(SDL_Event event, DrawingWindow &window) {
	if (event.type == SDL_KEYDOWN) {
		if (event.key.keysym.sym == SDLK_LEFT) std::cout << "LEFT" << std::endl;
		else if (event.key.keysym.sym == SDLK_RIGHT) std::cout << "RIGHT" << std::endl;
		else if (event.key.keysym.sym == SDLK_UP) std::cout << "UP" << std::endl;
		else if (event.key.keysym.sym == SDLK_DOWN) std::cout << "DOWN" << std::endl;
	} else if (event.type == SDL_MOUSEBUTTONDOWN) {
		window.savePPM("output.ppm");
		window.saveBMP("output.bmp");
	}
}

int main(int argc, char *argv[]) {
	DrawingWindow window = DrawingWindow(WIDTH, HEIGHT, false);
	SDL_Event event;
	
	float vertexScale = 0.17;

	std::unordered_map<std::string, Colour> materials = loadMtlFile("cornell-box.mtl");

	std::vector<ModelTriangle> faces = loadObjFile("cornell-box.obj", vertexScale, materials);
	glm::vec3 cameraPosition = glm::vec3(0.0, 0.0, 4.0);
	float focalLength = 2.0;

	while (true) {
		if (window.pollForInputEvents(event)){
			handleEvent(event, window);
			pointcloud(window, cameraPosition, faces, focalLength);
		}
		window.renderFrame();
	}
}