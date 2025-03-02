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
// directional light
Vec3f direct_light = Vec3f(0, 0, -1);
const int width  = 800;
const int height = 800;

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
void triangle(int i, float *zbuffer, TGAImage &image) {
    /*
     * rasterize i-th face from mesh
     */
    // fill in i-th face information, including pts, uvs
    Vec3f screen_coords[3];
    Vec3f world_coords[3];
    Vec2f uv_coords[3];
    for (int j = 0; j < 3; ++j) {
        world_coords[j] = model->vert(i, j);
        screen_coords[j] = world2screen(world_coords[j], width, height);
        uv_coords[j] = model->uv(i, j);
    }
    // CCW -> outside direction: normal = AB ^ AC
    Vec3f normal = (world_coords[1] - world_coords[0]) ^ (world_coords[2] - world_coords[0]);
    normal.normalize();
    float cos_theta = normal * (direct_light * -1.f);
    // backculling return directly
    if (cos_theta <= 0) {
        return;
    }
    // find bounding box of a triangle
    Vec2f bboxmin(image.width() - 1, image.height() - 1);
    Vec2f bboxmax(0, 0);
    for (int i = 0; i < 3; ++i) {
        bboxmin.x = std::max(0.f, std::min(bboxmin.x, screen_coords[i].x));
        bboxmin.y = std::max(0.f, std::min(bboxmin.y, screen_coords[i].y));
        bboxmax.x = std::min(image.width() - 1.f, std::max(bboxmax.x, screen_coords[i].x));
        bboxmax.y = std::min(image.height() - 1.f, std::max(bboxmax.y, screen_coords[i].y));
    }
    Vec3f P;
    for (P.y = bboxmin.y; P.y <= bboxmax.y; P.y++) {
        for (P.x = bboxmin.x; P.x <= bboxmax.x ; P.x++) {
            Vec3f bary_coords = barycentric(screen_coords, P);
            if (bary_coords.x < 0. || bary_coords.y < 0. || bary_coords.z < 0.) continue;
            // depth interpolation
            P.z = bary_coords.x * screen_coords[0].z + bary_coords.y * screen_coords[1].z + bary_coords.z * screen_coords[2].z;
            // pytorch3d frame like
            if (P.z > zbuffer[int(P.y * width + P.x)]) {
                zbuffer[int(P.y * width + P.x)] = P.z;
                Vec2f uv_P = uv_coords[0] * bary_coords.x + uv_coords[1] * bary_coords.y + uv_coords[2] * bary_coords.z;
                image.set(P.x, P.y, model->get_color(uv_P) * cos_theta);
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

    TGAImage image(width, height, TGAImage::RGB);
    // zbuffer for min depth recording
    float *zbuffer = new float[width * height];
    // initialization
    for (int i = 0; i < width * height; ++i) {
//        zbuffer[i] = -std::numeric_limits<float>::max();
        zbuffer[i] = std::numeric_limits<float>::lowest();
    }

    for (int i = 0; i < model->nfaces(); ++i) {
        triangle(i, zbuffer, image);
    }

    image.write_tga_file("rasterization_w_texture.tga");
    delete model;
    return 0;
}

