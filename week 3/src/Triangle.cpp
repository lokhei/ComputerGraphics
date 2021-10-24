#include <CanvasTriangle.h>
#include <DrawingWindow.h>
#include <Utils.h>
#include <fstream>
#include <vector>
#include <CanvasPoint.h>
#include <Colour.h>
#include <CanvasTriangle.h>
#include <TextureMap.h>
#include <TexturePoint.h>

#define WIDTH 320
#define HEIGHT 240

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
	for (int i=0; i<numberOfValues; i++) {
		points.push_back(CanvasPoint(round(xs[i]), round(ys[i])));
	}
	return points;
}


void drawLine(DrawingWindow &window, CanvasPoint from, CanvasPoint to, Colour colour) {
	float numberOfSteps = std::max(abs(to.x - from.x), abs(to.y - from.y));
	std :: vector<CanvasPoint> points = interpolateRoundPoints(from, to, numberOfSteps + 1);
	for (int i=0; i<numberOfSteps; i++) {
		window.setPixelColour(points[i].x, points[i].y, (255 << 24) + (colour.red << 16) + (colour.green << 8) + colour.blue);
	}
}


void drawStrokedTriangle(DrawingWindow &window, CanvasTriangle triangle, Colour colour) {
	drawLine(window, triangle[0], triangle[1], colour);
	drawLine(window, triangle[0], triangle[2], colour);
	drawLine(window, triangle[1], triangle[2], colour);
}


CanvasTriangle randomTriangle(DrawingWindow &window){
	CanvasTriangle triangle = CanvasTriangle(
		CanvasPoint(rand()%window.width, rand()%window.height),
		CanvasPoint(rand()%window.width, rand()%window.height),
		CanvasPoint(rand()%window.width, rand()%window.height)
	);
	return triangle;
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


void drawFilledTriangle(DrawingWindow &window, CanvasTriangle triangle, Colour colour) {
	sortTriangle(triangle);
	std :: vector<CanvasPoint> start = interpolateRoundPoints(triangle[0], triangle[1], triangle[1].y - triangle[0].y + 1);
	start.pop_back(); // Last row duplicated by next triangle so pop
	std :: vector<CanvasPoint> start2 = interpolateRoundPoints(triangle[1], triangle[2], triangle[2].y - triangle[1].y + 1);
	start.insert(start.end(), start2.begin(), start2.end());
	std :: vector<CanvasPoint> end = interpolateRoundPoints(triangle[0], triangle[2], triangle[2].y - triangle[0].y + 1);

	// Draw the filled in triangle
	for (int i=0; i<=triangle[2].y - triangle[0].y; i++) {
		drawLine(window, start[i], end[i], colour);
	}

	drawStrokedTriangle(window, triangle, Colour(255,255,255)); //outline

}


void drawTexturedTriangle(DrawingWindow &window, CanvasTriangle triangle, TextureMap &texMap) {
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

		for (float j=0.0; j<=numberOfSteps; j++) {
			uint32_t col = texMap.pixels[(round(texPoints[j].y) * texMap.width) + round(texPoints[j].x)];
			window.setPixelColour(points[j].x, points[j].y, col);
		}

	}
	
	
	drawStrokedTriangle(window, triangle, Colour(255, 255, 255));
}


void handleEvent(SDL_Event event, DrawingWindow &window) {
	if (event.type == SDL_KEYDOWN) {
		if (event.key.keysym.sym == SDLK_LEFT) std::cout << "LEFT" << std::endl;
		else if (event.key.keysym.sym == SDLK_RIGHT) std::cout << "RIGHT" << std::endl;
		else if (event.key.keysym.sym == SDLK_UP) std::cout << "UP" << std::endl;
		else if (event.key.keysym.sym == SDLK_DOWN) std::cout << "DOWN" << std::endl;
		else if (event.key.keysym.sym == SDLK_u) drawStrokedTriangle(window, randomTriangle(window), Colour(rand()%256, rand()%256, rand()%256));
		else if (event.key.keysym.sym == SDLK_f) drawFilledTriangle(window, randomTriangle(window), Colour(rand()%256, rand()%256, rand()%256));
	} else if (event.type == SDL_MOUSEBUTTONDOWN) {
		window.savePPM("output.ppm");
		window.saveBMP("output.bmp");
	}
}

int main(int argc, char *argv[]) {
	DrawingWindow window = DrawingWindow(WIDTH, HEIGHT, false);
	SDL_Event event;
	TextureMap textureMap = TextureMap("texture.ppm");
	CanvasTriangle triangle = CanvasTriangle(CanvasPoint(160, 10), CanvasPoint(300, 230), CanvasPoint(10, 150));
	triangle[0].texturePoint = TexturePoint(195, 5);
	triangle[1].texturePoint = TexturePoint(395, 380);
	triangle[2].texturePoint = TexturePoint(65, 330);
	while (true) {
		if (window.pollForInputEvents(event)) handleEvent(event, window);
		drawTexturedTriangle(window, triangle, textureMap);
		window.renderFrame();
	}

}
