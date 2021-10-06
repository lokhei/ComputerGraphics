#include <CanvasTriangle.h>
#include <DrawingWindow.h>
#include <Utils.h>
#include <fstream>
#include <vector>
#include <glm/glm.hpp>

#define WIDTH 320
#define HEIGHT 240


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

std :: vector<float> interpolateSingleFloats(float from, float to, int numberOfValues) {
	float step = (to - from) / (numberOfValues -1);
	std :: vector<float> v;
	for (int i = 0; i < numberOfValues; i++){
		v.push_back(from + i*step);
	}
	return v;
}

std :: vector<glm::vec3> interpolateThreeElementValues(glm::vec3 from, glm::vec3 to, int numberOfValues){
	int divisor = numberOfValues -1;
	glm::vec3 step = (to - from) / glm::vec3(divisor);
	std :: vector<glm::vec3> v;
	for (int i = 0; i < numberOfValues; i++){
		v.push_back(from + glm::vec3(i)*step);
	}
	return v;

}

void drawGreyScale(DrawingWindow &window){
	window.clearPixels();
	std :: vector<float> grey = interpolateSingleFloats(255, 0, window.width);
	for (size_t x = 0; x < window.width; x++) {
		float greyCol = grey[x];
		for (size_t y = 0; y < window.height; y++) {
			uint32_t colour = (255 << 24) + (int(greyCol) << 16) + (int(greyCol) << 8) + int(greyCol);
			window.setPixelColour(x, y, colour);
		}
	}
}

void drawColourScale(DrawingWindow &window){
	window.clearPixels();
	glm::vec3 topLeft(255, 0, 0);        // red 
	glm::vec3 topRight(0, 0, 255);       // blue 
	glm::vec3 bottomRight(0, 255, 0);    // green 
	glm::vec3 bottomLeft(255, 255, 0);   // yellow

	std :: vector<glm::vec3> leftColumn = interpolateThreeElementValues(topLeft, bottomLeft, window.height);
	std :: vector<glm::vec3> rightColumn = interpolateThreeElementValues(topRight, bottomRight, window.height);

	for (size_t y = 0; y < window.height; y++) {
		std :: vector<glm::vec3> pixelRow = interpolateThreeElementValues(leftColumn[y], rightColumn[y], window.width);
		for (size_t x = 0; x < window.width; x++) {
			uint32_t colour = (255 << 24) + (int(pixelRow[x][0]) << 16) + (int(pixelRow[x][1]) << 8) + int(pixelRow[x][2]);
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
	
	/* test interpolateSingleFloats
	std::vector<float> result;
	result = interpolateSingleFloats(2.2, 8.5, 7);
	for(size_t i=0; i<result.size(); i++) std::cout << result[i] << " ";
	*/


	/* test interpolateThreeElementValues
	std::vector<glm::vec3> result;
	result = interpolateThreeElementValues(glm::vec3(1, 4, 9.2), glm::vec3(4, 1, 9.8), 4);
	for(size_t i=0; i<result.size(); i++){
		for (int j = 0; j < 3; j++) std::cout << result[i][j] << " ";
		std::cout << std::endl;
	}
	
	*/
	SDL_Event event;
	while (true) {
		// We MUST poll for events - otherwise the window will freeze !
		if (window.pollForInputEvents(event)) handleEvent(event, window);
		//draw(window);
		drawGreyScale(window);
		drawColourScale(window);
		// Need to render the frame at the end, or nothing actually gets shown on the screen !
		window.renderFrame();
	}

	
}
