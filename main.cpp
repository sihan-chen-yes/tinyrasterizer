//
// Created by csh on 2/26/25.
// tinyrasterizer
//

#include <vector>
#include <cmath>
#include "tgaimage.h"
#include "model.h"
#include "geometry.h"

constexpr TGAColor white   = {255, 255, 255, 255}; // attention, BGRA order
constexpr TGAColor green   = {  0, 255,   0, 255};
constexpr TGAColor red     = {  0,   0, 255, 255};
constexpr TGAColor blue    = {255, 128,  64, 255};
constexpr TGAColor yellow  = {  0, 200, 255, 255};
Model *model = NULL;
const int width  = 800;
const int height = 800;

// naive method need to specify step size
//void line(int x0, int y0, int x1, int y1, TGAImage &image, TGAColor color) {
//    for (float t = 0.0; t < 1.0; t += 0.01) {
//        // linear interpolation
//        int x = x0 + t * (x1 - x0);
//        int y = y0 + t * (y1 - y0);
//        image.set(x, y, color);
//    }
//}

// Bresenham Algorithm
//void line(int x0, int y0, int x1, int y1, TGAImage &image, TGAColor color) {
//    bool steep = false;
//    // assuming x direction iteration
//    if (std::abs(x0 - x1) < std::abs(y0 - y1)) {
//        // transpose picture
//        std::swap(x0, y0);
//        std::swap(x1, y1);
//        steep = true;
//    }
//    // assuming x0 < x1
//    if (x0 > x1) {
//        std::swap(x0, x1);
//        std::swap(y0, y1);
//    }
//    for (int x = x0; x <= x1; ++x) {
//        float t = (float) (x - x0) / (x1 - x0);
//        int y = y0 + t * (y1 - y0);
//        if (!steep) {
//            // didn't transpose
//            image.set(x, y, color);
//        } else{
//            // remember transpose back
//            image.set(y, x, color);
//        }
//    }
//}

// Bresenham Algorithm Optimization v1
// convert division to add and multiply
//void line(int x0, int y0, int x1, int y1, TGAImage &image, TGAColor color) {
//    bool steep = false;
//    // assuming x direction iteration
//    if (std::abs(x0 - x1) < std::abs(y0 - y1)) {
//        // transpose picture
//        std::swap(x0, y0);
//        std::swap(x1, y1);
//        steep = true;
//    }
//    // assuming x0 < x1
//    if (x0 > x1) {
//        std::swap(x0, x1);
//        std::swap(y0, y1);
//    }
//    // d_x = 1, d_y = k x d_x = k
//    int dy = y1 - y0;
//    int dx = x1 - x0;
//    float d_y_step = std::abs((float) dy / dx);
//    float d_y_total = 0;
//    int y = y0;
//    for (int x = x0; x <= x1; ++x) {
//        if (!steep) {
//            // didn't transpose
//            image.set(x, y, color);
//        } else{
//            // remember transpose back
//            image.set(y, x, color);
//        }
//        d_y_total += d_y_step;
//        // beyond half grid
//        if (d_y_total > 0.5) {
//            y += y1 > y0 ? 1 : -1;
//            // should move whole grid (0.5 * 2)
//            d_y_total -= 1;
//        }
//    }
//}

// Bresenham Algorithm Optimization v2
// remove float completely via transforming the inequality (multiply 2 on both sides)
void line(int x0, int y0, int x1, int y1, TGAImage &image, TGAColor color) {
    bool steep = false;
    // assuming x direction iteration
    if (std::abs(x0 - x1) < std::abs(y0 - y1)) {
        // transpose picture
        std::swap(x0, y0);
        std::swap(x1, y1);
        steep = true;
    }
    // assuming x0 < x1
    if (x0 > x1) {
        std::swap(x0, x1);
        std::swap(y0, y1);
    }
    // d_x = 1, d_y = k x d_x = k
    int dy = y1 - y0;
    int dx = x1 - x0;
    int d_error = 2 * std::abs(dy);
    int d_error_total = 0;
    int y = y0;
    for (int x = x0; x <= x1; ++x) {
        if (!steep) {
            // didn't transpose
            image.set(x, y, color);
        } else{
            // remember transpose back
            image.set(y, x, color);
        }
        d_error_total += d_error;
        // beyond half grid
        if (d_error_total > dx) {
            y += y1 > y0 ? 1 : -1;
            // should move whole grid (dx * 2)
            d_error_total -= 2 * dx;
        }
    }
}

int main(int argc, char** argv) {
    if (2==argc) {
        model = new Model(argv[1]);
    } else {
        model = new Model("obj/african_head.obj");
    }

    TGAImage image(width, height, TGAImage::RGB);
    for (int i=0; i<model->nfaces(); i++) {
        std::vector<int> face = model->face(i);
        for (int j=0; j<3; j++) {
            Vec3f v0 = model->vert(face[j]);
            Vec3f v1 = model->vert(face[(j+1)%3]);
            // transform from NDC to img coordinates
            int x0 = (v0.x+1.)*width/2.;
            int y0 = (v0.y+1.)*height/2.;
            int x1 = (v1.x+1.)*width/2.;
            int y1 = (v1.y+1.)*height/2.;
            line(x0, y0, x1, y1, image, white);
        }
    }

    image.write_tga_file("output.tga");
    delete model;
    return 0;
}