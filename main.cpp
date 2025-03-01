//
// Created by csh on 2/26/25.
// tinyrasterizer
//

#include <vector>
#include <cmath>
#include <iostream>
#include "tgaimage.h"
#include "model.h"
#include "geometry.h"

TGAColor white   = TGAColor (255, 255, 255, 255); // attention, BGRA order
TGAColor green   = TGAColor (  0, 255,   0, 255);
TGAColor red     = TGAColor (  255,   0, 0, 255);
TGAColor blue    = TGAColor (64, 128,  255, 255);
TGAColor yellow  = TGAColor (  255, 200, 0, 255);
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
void line(Vec2i p0, Vec2i p1, TGAImage &image, TGAColor color) {
    int x0 = p0.x;
    int y0 = p0.y;
    int x1 = p1.x;
    int y1 = p1.y;
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

// Wireframe rendering
//int main(int argc, char** argv) {
//    if (2 == argc) {
//        model = new Model(argv[1]);
//    } else {
//        model = new Model("obj/african_head.obj");
//    }
//
//    TGAImage image(width, height, TGAImage::RGB);
//    for (int i =0 ; i < model->nfaces(); i++) {
//        std::vector<int> face = model->face(i);
//        for (int j = 0; j < 3; j++) {
//            Vec3f v0 = model->vert(face[j]);
//            Vec3f v1 = model->vert(face[(j + 1) % 3]);
//            // transform from NDC to img coordinates
//            line(Vec2i((v0.x + 1.) * width / 2., (v0.y + 1.) * height / 2.),
//                 Vec2i((v1.x + 1.) * width / 2., (v1.y + 1.) * height / 2.), image, white);
//        }
//    }
//
//    image.write_tga_file("Bresenham.tga");
//    delete model;
//    return 0;
//}

// Naive sweeping line Algorithm
//void triangle(Vec2i t0, Vec2i t1, Vec2i t2, TGAImage &image, TGAColor color) {
//    // sort according y coordinate
//    if (t0.y > t1.y) std::swap(t0, t1);
//    if (t0.y > t2.y) std::swap(t0, t2);
//    if (t1.y > t2.y) std::swap(t1, t2);
//    float total_height = t2.y - t0.y;
//    float segment_height = t1.y - t0.y + 1;
//    // the first half of triangle
//    for (int y = t0.y; y <= t1.y; ++y) {
//        float alpha = (float) (y - t0.y) / total_height;
//        float beta = (float) (y - t0.y) / segment_height;
//        Vec2i A = t0 + (t2 - t0) * alpha;
//        Vec2i B = t0 + (t1 - t0) * beta;
//        // assuming A is on the left side and B is on the right side
//        if (A.x > B.x) std::swap(A, B);
//        // sweeping each row
//        for (int x = A.x; x <= B.x; ++x) {
//            image.set(x, y, color);
//        }
//    }
//    // the second half of triangle
//    segment_height = t2.y - t1.y + 1;
//    for (int y = t1.y; y <= t2.y; ++y) {
//        float alpha = (float) (y - t0.y) / total_height;
//        float beta = (float) (y - t1.y) / segment_height;
//        Vec2i A = t0 + (t2 - t0) * alpha;
//        Vec2i B = t1 + (t2 - t1) * beta;
//        // assuming A is on the left side and B is on the right side
//        if (A.x > B.x) std::swap(A, B);
//        // sweeping each row
//        for (int x = A.x; x <= B.x; ++x) {
//            image.set(x, y, color);
//        }
//    }
//}

Vec3f barycentric(Vec2i *pts, Vec2i P) {
    // *pts: A B C
    // [u, v, 1] will be orthogonal to vec1 and vec2
    // [AB,AC,PA].x X [AB,AC,PA].y
    Vec3f n = Vec3f(pts[1].x - pts[0].x, pts[2].x - pts[0].x, pts[0].x - P.x) ^
              Vec3f(pts[1].y - pts[0].y, pts[2].y - pts[0].y, pts[0].y - P.y);
    // if degrade, return fix invalid value
    if (fabs(n.z) < EPSILON) return Vec3f(-1., 1., 1.);
    // normalize and return
    return Vec3f(1. - (n.x + n.y) / n.z, n.x / n.z, n.y / n.z);
}

// drawing triangles iterating bbox
void triangle(Vec2i *pts, TGAImage &image, TGAColor color) {
    Vec2i bboxmin(image.height() - 1, image.width() - 1);
    Vec2i bboxmax(0, 0);
    for (int i = 0; i < 3; ++i) {
        bboxmin.x = std::max(0, std::min(bboxmin.x, pts[i].x));
        bboxmin.y = std::max(0, std::min(bboxmin.y, pts[i].y));
        bboxmax.x = std::min(image.height() - 1, std::max(bboxmax.x, pts[i].x));
        bboxmax.y = std::min(image.width() - 1, std::max(bboxmax.y, pts[i].y));
    }
    Vec2i P;
    for (P.y = bboxmin.y; P.y <= bboxmax.y; P.y++) {
        for (P.x = bboxmin.x; P.x <= bboxmax.x ; P.x++) {
            Vec3f bary_coords = barycentric(pts, P);
            if (bary_coords.x < 0. || bary_coords.y < 0. || bary_coords.z < 0.) continue;
            image.set(P.x, P.y, color);
        }
    }
}

//int main(int argc, char** argv) {
//    if (2 == argc) {
//        model = new Model(argv[1]);
//    } else {
//        model = new Model("obj/african_head.obj");
//    }
//
//    TGAImage image(width, height, TGAImage::RGB);
//    for (int i = 0; i < model->nfaces(); ++i) {
//        std::vector<int> face = model->face(i);
//        Vec2i screen_coords[3];
//        for (int j = 0; j < 3; ++j) {
//            Vec3f world_coords = model->vert(face[j]);
//            screen_coords[j] = Vec2i((world_coords.x + 1) * image.width() / 2.,
//                                     (world_coords.y + 1) * image.height() / 2.);
//        }
//        TGAColor rand_color = TGAColor(rand() % 255, rand() % 255, rand() % 255, 255);
//        triangle(screen_coords, image, rand_color);
//    }
//
//    image.write_tga_file("triangle.tga");
//    delete model;
//    return 0;
//}


int main(int argc, char** argv) {
    if (2 == argc) {
        model = new Model(argv[1]);
    } else {
        model = new Model("obj/african_head.obj");
    }

    Vec3f light_dir(0,0,1); // define light direction

    TGAImage image(width, height, TGAImage::RGB);
    for (int i = 0; i < model->nfaces(); ++i) {
        std::vector<int> face = model->face(i);
        Vec2i screen_coords[3];
        Vec3f world_coords[3];
        for (int j = 0; j < 3; ++j) {
            world_coords[j] = model->vert(face[j]);
            screen_coords[j] = Vec2i((world_coords[j].x + 1) * image.width() / 2.,
                                     (world_coords[j].y + 1) * image.height() / 2.);
        }

        // compute face normal
        Vec3f n = (world_coords[2] - world_coords[0]) ^ (world_coords[1] - world_coords[0]);
        n.normalize();
        // reverse light direction and calculate cos theta
        float cosTheta = n * (light_dir * -1);
        if (cosTheta > 0) {
            triangle(screen_coords, image, TGAColor(cosTheta * 255, cosTheta * 255, cosTheta * 255, 255));
        }
    }

    image.write_tga_file("backculling.tga");
    delete model;
    return 0;
}