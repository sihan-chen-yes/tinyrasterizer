//
// Created by csh on 3/15/25.
//

#ifndef TINYRASTERIZER_OUR_GL_H
#define TINYRASTERIZER_OUR_GL_H

#include "tgaimage.h"
#include "geometry.h"

extern Matrix ModelView;
extern Matrix Viewport;
extern Matrix Projection;
const float max_depth = 2000.f;

void lookat(Vec3f eye, Vec3f center, Vec3f up);
void projection(float coeff=0.f); // coeff = -1/c
void viewport(int x, int y, int w, int h);

struct IShader {
    virtual ~IShader();
    virtual Vec4f vertex(int iface, int jvert) = 0;
    virtual bool fragment(Vec3f bar, TGAColor &color) = 0;
};

void triangle(Vec4f *h_pts, IShader &shader, TGAImage &image, float *zbuffer);
#endif //TINYRASTERIZER_OUR_GL_H
