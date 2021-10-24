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

#define WIDTH 640
#define HEIGHT 480

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
		} else if(line[0] == "map_Kd") {
			Colour colour = colours[colour_name];
			colour.name = line[1];
			colours[colour_name] = colour;
		}
	}
	inputStream.close();
	return colours;
}

std::vector<ModelTriangle> loadObjFile(const std::string &filename, float scale) {
	std::vector<glm::vec3> vertices;
	std::vector<ModelTriangle> faces;
	std::vector<TexturePoint> texturePoints;
	std::unordered_map<std::string, Colour> materials;

	std::ifstream inputStr(filename, std::ifstream::in);
	std::string nextLine;
	Colour colour;
	while (std::getline(inputStr, nextLine)) { //extracts from inputStr and stores into nextLine
		std::vector<std::string> vector = split(nextLine, ' '); //split line by spaces
		if (vector[0] == "mtllib") {
			materials = loadMtlFile(vector[1]);
		}else if (vector[0] == "usemtl") {
			colour =materials[vector[1]];
		}else if (vector[0] == "v") {
			vertices.push_back(glm::vec3(
				std::stof(vector[1]) * scale, //string to float
				std::stof(vector[2]) * scale,
				std::stof(vector[3]) * scale
			));
		}else if (vector[0] == "f") { //indexed from 1
			std::vector<std::string> l1 = split(vector[1], '/');
			std::vector<std::string> l2 = split(vector[2], '/');
			std::vector<std::string> l3 = split(vector[3], '/');
			ModelTriangle triangle(
				vertices[std::stoi(l1[0])-1], 
				vertices[std::stoi(l2[0])-1], 
				vertices[std::stoi(l3[0])-1], 
				colour);
			if(l1[1] != "") {
				triangle.texturePoints[0] = texturePoints[stoi(l1[1])-1];
				triangle.texturePoints[1] = texturePoints[stoi(l2[1])-1];
				triangle.texturePoints[2] = texturePoints[stoi(l3[1])-1];
			} 
			faces.push_back(triangle);
		}else if(vector[0] == "vt") {
			texturePoints.push_back(TexturePoint(stof(vector[1]), stof(vector[2])));
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

std::vector<CanvasPoint> interpolateRoundPoints(CanvasPoint from, CanvasPoint to, int numberOfValues) {
	std::vector<CanvasPoint> interpolatedValues;
	auto xs = interpolateSingleFloats(from.x, to.x, numberOfValues);
	auto ys = interpolateSingleFloats(from.y, to.y, numberOfValues);
	auto depths = interpolateSingleFloats(from.depth, to.depth, numberOfValues);
	for (int i=0; i<numberOfValues; i++) {
		interpolatedValues.push_back(CanvasPoint(round(xs[i]), round(ys[i]), depths[i]));
	}
	return interpolatedValues;
}

std::vector<CanvasPoint> interpolateRoundPoints(TexturePoint from, TexturePoint to, int numberOfValues) {
	std::vector<CanvasPoint> interpolatedValues;
	auto xs = interpolateSingleFloats(from.x, to.x, numberOfValues);
	auto ys = interpolateSingleFloats(from.y, to.y, numberOfValues);
	for (int i=0; i<numberOfValues; i++) {
		interpolatedValues.push_back(CanvasPoint(round(xs[i]), round(ys[i])));
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
	float numberOfSteps = std::max(std::max(abs(to.x - from.x), abs(to.y - from.y)), 1.0f);
	std :: vector<CanvasPoint> points = interpolateRoundPoints(from, to, numberOfSteps + 1);
	for (int i=0; i<=numberOfSteps; i++) {
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

void drawFilledTriangle(DrawingWindow &window, CanvasTriangle triangle, Colour colour, std::vector<std::vector<float>> &depthBuffer, bool outline) {
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


void drawTexturedTriangle(DrawingWindow &window, CanvasTriangle triangle, TextureMap &texMap, std::vector<std::vector<float>> &depthBuffer, bool outline) {
	sortTriangle(triangle);
	std :: vector<CanvasPoint> xStart = interpolateRoundPoints(triangle[0], triangle[1], triangle[1].y - triangle[0].y + 1);
	xStart.pop_back();
	std :: vector<CanvasPoint> xStart2 = interpolateRoundPoints(triangle[1], triangle[2], triangle[2].y - triangle[1].y + 1);
	xStart.insert(xStart.end(), xStart2.begin(), xStart2.end());
	std :: vector<CanvasPoint> xEnd = interpolateRoundPoints(triangle[0], triangle[2], triangle[2].y - triangle[0].y + 1);


	std :: vector<CanvasPoint> textureStarts = interpolateRoundPoints(triangle[0].texturePoint, triangle[1].texturePoint, triangle[1].y - triangle[0].y + 1);
	textureStarts.pop_back();
	std :: vector<CanvasPoint> textureStarts2 = interpolateRoundPoints(triangle[1].texturePoint, triangle[2].texturePoint, triangle[2].y - triangle[1].y + 1);
	textureStarts.insert(textureStarts.end(), textureStarts2.begin(), textureStarts2.end());
	std :: vector<CanvasPoint> textureEnds = interpolateRoundPoints(triangle[0].texturePoint, triangle[2].texturePoint, triangle[2].y - triangle[0].y + 1);


	for(int i=0; i <= triangle[2].y - triangle[0].y; i++) {
		float numberOfSteps = abs(xStart[i].x - xEnd[i].x);
		//interpolate between start and end of rake to find texture point for pixel
		std :: vector<CanvasPoint> texPoints = interpolateRoundPoints(textureStarts[i], textureEnds[i], numberOfSteps + 1);
		std :: vector<CanvasPoint> points = interpolateRoundPoints(xStart[i], xEnd[i], numberOfSteps + 1);

		for (float j=0.0; j <= numberOfSteps; j++) {
			int x = points[j].x;
			int y = points[j].y;
			if(x >= 0 && x < window.width && y >= 0 && y < window.height) {
				float pointDepth = 1 / -points[i].depth;	
				if (pointDepth > depthBuffer[y][x]) {
					depthBuffer[y][x] = pointDepth;
					uint32_t col = texMap.pixels[texPoints[j].y * texMap.width + texPoints[j].x];
					window.setPixelColour(points[j].x, points[j].y, col);
				}
			}
		}
	}

	if (outline) drawStrokedTriangle(window, triangle, Colour(255, 255, 255), depthBuffer);
}



CanvasPoint getCanvasIntersectionPoint(DrawingWindow &window, glm::vec3 cameraPosition, glm::mat3 camOrientation, glm::vec3 vertexPosition, float focalLength){

	int planeMultiplier = 600;
	glm::vec3 vertex = (vertexPosition - cameraPosition) * camOrientation;
	float u = -round(planeMultiplier*focalLength * (vertex.x / vertex.z)) + (window.width / 2);
	float v = round(planeMultiplier*focalLength * (vertex.y / vertex.z)) + (window.height / 2);
	float z = vertex.z;

	return CanvasPoint(u,v,z);
}

void drawObj(DrawingWindow &window, glm::vec3 cameraPosition, glm::mat3 camOrientation, std::vector<ModelTriangle> faces, float focalLength, TextureMap textureMap) {
	window.clearPixels();
	std::vector<std::vector<float>> depthBuffer (window.height, std::vector<float>(window.width, 0));

	for (int i=0; i<faces.size(); i++) {
		ModelTriangle face = faces[i];
		CanvasTriangle triangle = CanvasTriangle();
		for (int j=0; j<face.vertices.size(); j++) {
			glm::vec3 vertexPosition = face.vertices[j];
			triangle.vertices[j] = getCanvasIntersectionPoint(window,cameraPosition,camOrientation, vertexPosition, focalLength );
			triangle.vertices[j].texturePoint = face.texturePoints[j];
		}

		if (face.colour.name != "") {

			TextureMap texture(face.colour.name);
			for(int j = 0; j < triangle.vertices.size(); j++) {
				triangle.vertices[j].texturePoint.x = fmod(triangle.vertices[j].texturePoint.x, 1 ) * texture.width;
				triangle.vertices[j].texturePoint.y = texture.height-fmod(triangle.vertices[j].texturePoint.y, 1) * texture.height;
			}

			drawTexturedTriangle(window, triangle, textureMap, depthBuffer, false);
		}
		else {
			drawFilledTriangle(window, triangle, faces[i].colour, depthBuffer, false);
		}
	}
	
}

template <typename T>
void rotateY(T &camera, float angle) {
	glm::mat3 rotationMatrix = glm::mat3(
		cos(angle), 0.0, -sin(angle),
		0.0, 1.0, 0.0,
		sin(angle), 0.0, cos(angle)
	);
	camera = rotationMatrix * camera;
}

template <typename T>
void rotateX(T &camera, float angle) {
	glm::mat3 rotationMatrix = glm::mat3(
		1.0, 0.0, 0.0,
		0.0, cos(angle), sin(angle),
		0.0, -sin(angle), cos(angle)
	);
	camera = rotationMatrix * camera;
}


void handleEvent(SDL_Event event, DrawingWindow &window, glm::vec3 &camPos, glm::mat3 &camOri) {
	if (event.type == SDL_KEYDOWN) {
		//camera movement
		if (event.key.keysym.sym == SDLK_e) camPos.y += 0.1; //up
		else if (event.key.keysym.sym == SDLK_q) camPos.y -= 0.1; //down
		else if (event.key.keysym.sym == SDLK_w) camPos.z -= 0.1; //forwards
		else if (event.key.keysym.sym == SDLK_s) camPos.z += 0.1; //backwards
		else if (event.key.keysym.sym == SDLK_d) camPos.x += 0.1; //right
		else if (event.key.keysym.sym == SDLK_a) camPos.x -= 0.1; //left
		//rotation
		else if (event.key.keysym.sym == SDLK_UP) rotateY(camPos, 0.1); //rotate Y C
		else if (event.key.keysym.sym == SDLK_RIGHT) rotateX(camPos, 0.1); //rotate X C
		else if (event.key.keysym.sym == SDLK_DOWN) rotateY(camPos, -0.1); //rotate Y antiC
		else if (event.key.keysym.sym == SDLK_LEFT) rotateX(camPos, -0.1); //rotate X antiC
		//orientation
		else if (event.key.keysym.sym == SDLK_i) rotateY(camOri, 0.1); //panning
		else if (event.key.keysym.sym == SDLK_k) rotateY(camOri, -0.1);
		else if (event.key.keysym.sym == SDLK_l) rotateX(camOri, 0.1); //tilting
		else if (event.key.keysym.sym == SDLK_j) rotateX(camOri, -0.1);


	} else if (event.type == SDL_MOUSEBUTTONDOWN) {
		window.savePPM("output.ppm");
		window.saveBMP("output.bmp");
	}
}

int main(int argc, char *argv[]) {
	DrawingWindow window = DrawingWindow(WIDTH, HEIGHT, false);
	SDL_Event event;

	TextureMap textureMap = TextureMap("texture.ppm");
	
	float vertexScale = 0.17;
	
	std::vector<ModelTriangle> faces = loadObjFile("textured-cornell-box.obj", vertexScale);
	glm::vec3 cameraPosition = glm::vec3(0.0, 0.0, 4.0);
	float focalLength = 2.0;

	glm::mat3 camOrientation(
		1.0, 0.0, 0.0,
		0.0, 1.0, 0.0,
		0.0, 0.0, 1.0
		);

	while (true) {
		if (window.pollForInputEvents(event)){
			handleEvent(event, window, cameraPosition, camOrientation);
			drawObj(window, cameraPosition, camOrientation, faces, focalLength, textureMap);
		}
		window.renderFrame();
	}
}