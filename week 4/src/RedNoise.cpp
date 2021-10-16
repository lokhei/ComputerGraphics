#include <CanvasTriangle.h>
#include <DrawingWindow.h>
#include <Utils.h>
#include <fstream>
#include <vector>
#include <ModelTriangle.h>

#define WIDTH 320
#define HEIGHT 240


std::vector<Colour> loadMtlFile(const std::string &filename) {
	std::vector<Colour> colours;

	std::ifstream inputStream(filename, std::ifstream::in);
	std::string nextLine;
	std::string name;

	while (std::getline(inputStream, nextLine)) {
		std::vector<std::string> line = split(nextLine, ' ');

		if (line[0] == "newmtl") {
			name = line[1];
		} else if (line[0] == "Kd") {
			colours.push_back(Colour(
				name,
				(int)(std::stof(line[1]) * 255),
				(int)(std::stof(line[2]) * 255),
				(int)(std::stof(line[3]) * 255)
			));
		}
	}
	return colours;
}


std::vector<ModelTriangle> loadObjFile(const std::string &filename, float scale, std::vector<Colour> materials) {
	std::vector<glm::vec3> vertices;
	std::vector<ModelTriangle> faces;

	std::ifstream inputStr(filename, std::ifstream::in);
	std::string nextLine;
	Colour colour;
	while (std::getline(inputStr, nextLine)) { //extracts from inputStr and stores into nextLine
		std::vector<std::string> vector = split(nextLine, ' '); //split line by spaces
		
		if (vector[0] == "usemtl") {
			for (int i=0; i<materials.size(); i++) {
				if (materials[i].name == vector[1]) {
					colour = materials[i];
					break;
				}
			}
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
	return faces;
}


void draw(DrawingWindow &window) {
	window.clearPixels();
	for (size_t y = 0; y < window.height; y++) {
		for (size_t x = 0; x < window.width; x++) {
			float red = rand() % 256;
			float green = 0.0;
			float blue = 0.0;
			uint32_t colour = (255 << 24) + (int(red) << 16) + (int(green) << 8) + int(blue);
			window.setPixelColour(x, y, colour);
		}
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

	std::vector<Colour> materials = loadMtlFile("cornell-box.mtl");
	std::vector<ModelTriangle> faces = loadObjFile("cornell-box.obj", vertexScale, materials);
	for (int i=0; i<faces.size(); i++) {
		std::cout << faces[i] << std::endl;
	}

	while (true) {
		// We MUST poll for events - otherwise the window will freeze !
		if (window.pollForInputEvents(event)) handleEvent(event, window);
		draw(window);
		// Need to render the frame at the end, or nothing actually gets shown on the screen !
		window.renderFrame();
	}
}
