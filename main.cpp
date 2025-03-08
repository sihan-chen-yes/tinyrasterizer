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
const int depth = 255;
// camera location
Vec3f camera(0, 0, 3);

Vec3f m2v(Matrix m) {
    // matrix form to vector form
    return Vec3f(m[0][0]/m[3][0], m[1][0]/m[3][0], m[2][0]/m[3][0]);
}

Matrix v2m(Vec3f v) {
    // vector form to matrix form
    Matrix m(4, 1);
    m[0][0] = v.x;
    m[1][0] = v.y;
    m[2][0] = v.z;
    m[3][0] = 1.f;
    return m;
}

Matrix get_view_matrix(int x0, int y0, int w, int h, int d) {
    /*
     * construct wolrd2view transformation matrix
     */
    Matrix m = Matrix::identity(4);
    m[0][3] = x0 + w / 2.f;
    m[1][3] = y0 + h / 2.f;
    m[2][3] = d / 2.f;

    m[0][0] = w / 2.f;
    m[1][1] = h / 2.f;
    m[2][2] = d / 2.f;
    return m;
}

Matrix get_projection_matrix(float camera_z) {
    /*
     * get perspective get_projection_matrix matrix
     */
    Matrix m = Matrix::identity(4);
    m[3][2] = -1.f / camera_z;
    return m;
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
    // simple orthographic get_projection_matrix
//    return Vec3f(int((v.x + 1) * width / 2 + 0.5), int((v.y + 1) * height / 2 + 0.5), v.z);

    // simple perspective get_projection_matrix
    Matrix view = get_view_matrix(width / 8, height / 8, width * 3 / 4, height * 3 / 4, depth);
    Matrix projection = get_projection_matrix(camera.z);
    // perspective projection in 3D space then to project to screen
    return m2v(view * projection * v2m(v));
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
    // round screen_coords to int to meet assumption
    for (int j = 0; j < 3; ++j) {
        screen_coords[j].x = float(int(screen_coords[j].x + 0.5));
        screen_coords[j].y = float(int(screen_coords[j].y + 0.5));
        screen_coords[j].z = float(int(screen_coords[j].z + 0.5));
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
        zbuffer[i] = -std::numeric_limits<float>::max();
    }

    for (int i = 0; i < model->nfaces(); ++i) {
        triangle(i, zbuffer, image);
    }

    image.write_tga_file("rasterization_perspective.tga");

    TGAImage zbimage(width, height, TGAImage::GRAYSCALE);
    // normalize depth, search for max and min valid depth
    float z_min = std::numeric_limits<float>::max();
    float z_max = -std::numeric_limits<float>::max();
    for (int x = 0; x < width; ++x) {
        for (int y = 0; y < height; ++y) {
            float z = zbuffer[x + y * width];
            // ensure z is valid, then record
            if (z > 0.f) {
                z_min = std::min(z, z_min);
                z_max = std::max(z, z_max);
            }
        }
    }
    std::cout << "z_min:" <<  z_min << std::endl;
    std::cout << "z_max:" << z_max << std::endl;
    for (int x = 0; x < width; x++) {
        for (int y = 0; y < height; y++) {
            float z = zbuffer[x + y * width];
            float z_normalized = (z - z_min) / (z_max - z_min);
            zbimage.set(x, y, white * z_normalized);
        }
    }
    zbimage.write_tga_file("zbuffer_debug.tga");

    delete model;
    return 0;
}
