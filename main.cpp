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

Vec3f barycentric(Vec3f *pts, Vec3f P) {
    /*
     * calculate the barycentric coordinates given a triangle ABC and point P
     */
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

Vec3f world2screen(Vec3f v, int width, int height) {
    /*
     * projecting world coordinates to screen coordinates
     */
    // 0.5 for rounding
    return Vec3f(int((v.x + 1) * width / 2 + 0.5), int((v.y + 1) * height / 2 + 0.5), v.z);
}

// drawing triangles iterating bbox
void triangle(Vec3f *pts, float *zbuffer, TGAImage &image, TGAColor* vert_color) {
    /*
     * pts: in screen coordinates
     */
    // find bounding box of a triangle
    Vec2f bboxmin(image.width() - 1, image.height() - 1);
    Vec2f bboxmax(0, 0);
    for (int i = 0; i < 3; ++i) {
        bboxmin.x = std::max(0.f, std::min(bboxmin.x, pts[i].x));
        bboxmin.y = std::max(0.f, std::min(bboxmin.y, pts[i].y));
        bboxmax.x = std::min(image.width() - 1.f, std::max(bboxmax.x, pts[i].x));
        bboxmax.y = std::min(image.height() - 1.f, std::max(bboxmax.y, pts[i].y));
    }
    Vec3f P;
    for (P.y = bboxmin.y; P.y <= bboxmax.y; P.y++) {
        for (P.x = bboxmin.x; P.x <= bboxmax.x ; P.x++) {
            Vec3f bary_coords = barycentric(pts, P);
            if (bary_coords.x < 0. || bary_coords.y < 0. || bary_coords.z < 0.) continue;
            // depth interpolation
            P.z = bary_coords.x * pts[0].z + bary_coords.y * pts[1].z + bary_coords.z * pts[2].z;
            // pytorch3d frame like
            if (P.z > zbuffer[int(P.y * width + P.x)]) {
                zbuffer[int(P.y * width + P.x)] = P.z;
                TGAColor color = vert_color[0] * bary_coords.x + vert_color[1] * bary_coords.y + vert_color[2] * bary_coords.z;
                image.set(P.x, P.y, color);
            }
        }
    }
}

void rasterize(Vec2i p0, Vec2i p1, TGAImage &image, TGAColor color, int ybuffer[]) {
    if (p0.x>p1.x) {
        std::swap(p0, p1);
    }
    for (int x=p0.x; x<=p1.x; x++) {
        float t = (x-p0.x)/(float)(p1.x-p0.x);
        int y = p0.y*(1.-t) + p1.y*t;
        if (ybuffer[x]<y) {
            ybuffer[x] = y;
            image.set(x, 0, color);
        }
    }
}

int main(int argc, char** argv) {
    if (2 == argc) {
        model = new Model(argv[1]);
    } else {
        model = new Model("obj/african_head.obj");
    }
    // directional light
    Vec3f direct_light = Vec3f(0, 0, -1);

    TGAImage image(width, height, TGAImage::RGB);
    // zbuffer for min depth recording
    float *zbuffer = new float[width * height];
    // initialization
    for (int i = 0; i < width * height; ++i) {
//        zbuffer[i] = -std::numeric_limits<float>::max();
        zbuffer[i] = std::numeric_limits<float>::lowest();
    }

    for (int i = 0; i < model->nfaces(); ++i) {
        Vec3f screen_coords[3];
        Vec3f world_coords[3];
        for (int j = 0; j < 3; ++j) {
            world_coords[j] = model->vert(i, j);
            screen_coords[j] = world2screen(world_coords[j], width, height);
        }
        // CCW -> outside direction: normal = AB ^ AC
        Vec3f normal = (world_coords[1] - world_coords[0]) ^ (world_coords[2] - world_coords[0]);
        normal.normalize();
        float cos_theta = normal * (direct_light * -1.f);

        if (cos_theta > 0) {
            TGAColor colors_vert[3];
            for (int j = 0; j < 3; ++j) {
                // first interpolate uv then query
                colors_vert[j] = TGAColor(255, 255, 255, 255) * cos_theta;
            }
            triangle(screen_coords, zbuffer, image, colors_vert);
        }
    }

    image.write_tga_file("rasterization_wo_texture.tga");
    delete model;
    return 0;
}

