# ComputerGraphics

COMS30020 Computer Graphics Coursework

Building a C++ renderer using a variety of different approaches to 3D rendering

## Animation

## Features Implemented
- [x] OBJ geometry and material file loading
- [x] Wireframe rendering
- [x] Flat shaded triangle rasterising
- [x] Texture Mapping
- [x] Keyboard control of camera position
- [x] Changing camera orientation (using orientation matrix)
- [x] Camera LookAt and Orbit
- [x] Ambient lighting
- [x] Diffuse lighting (proximity and angle-of-incidence)
- [x] Specular lighting
- [x] Hard Shadow
- [x] Gouraud shading
- [x] Phong Shading
- [x] Realistic soft shadows
- [x] Reflective materials - Mirrors
- [x] Refractive materials (e.g. glass, water etc.)
- [x] Export of image files and creation of animated video

### Keyboard Controls

`0` `1` `2` to switch between wireframe, rasterising and raytracing respectively

`3` `4` `5` to switch between normal, gouraud and phong lighting respectively

`6` `7` `8` to switch proximity, incidence and specular lighting on/off respectively

`q` `e` `w` `a` `s` `d` to translate the camera position down, up, forwards, left, back and right respectively

`left` `right` arrow keys to rotate the camera position about the X axis anti-clockwise and clockwise respectively

`up` `down` arrow keys to rotate the camera position about the Y axis clockwise and anti-clockwise respectively

`l` `j` keys to rotate the camera in the X axis (tilting)  respectively

`i` `k` keys to rotate the camera in the Y axis (panning)  respectively

`o` to start orbiting (about the y axis) around origin

`r` to reset the camera view

`g` `h` to switch hard ad soft shadows on/off respectively

Clicking the mouse inside the SDL window to generate PPM and BMP files of the current content of the window 

## Usage

### Build and run

To build and run this project, use the included `Makefile`

```bash
$ make
```

or for optimised build rule

```bash
$ make speedy
```
