//
// Created by csh on 3/15/25.
//

#include <cmath>
#include <limits>
#include <cstdlib>
#include "our_gl.h"

// Model-View-Projection-Viewport(MVPV) matrix
Matrix ModelView;
Matrix Viewport;
Matrix Projection;

IShader::~IShader() {}

void lookat(Vec3f eye, Vec3f center, Vec3f up) {
    /*
     * eye: camera location
     * center: looking at
     * up: up direction of camera
     *
     * return: View matrix to transform world coords to camera coords
     */

    // find camera local basis to build rotation matrix
    Vec3f z = (eye - center).normalize();
    // up and z share the same plane
    Vec3f x = (up ^ z).normalize();
    // up and y not necessarily aligned
    Vec3f y = (z ^ x).normalize();
    // build transformation matrix
    // M_inv = M.T
    // R|T inv ==> R.T|-R.T x T
    Matrix M_inv = Matrix::identity(4);
    Matrix Tr = Matrix::identity(4);
    for (int i = 0; i < 3; ++i) {
        M_inv[0][i] = x[i];
        M_inv[1][i] = y[i];
        M_inv[2][i] = z[i];
        Tr[i][3] = -center[i];
    }
    ModelView = M_inv * Tr;
}

void projection(float coeff) {
    /*
     * construct Projection transformation matrix
     * coeff = -1 / c
     * c = (eye - center).norm()
     */
    Projection = Matrix::identity(4);
    Projection[3][2] = coeff;
}

void viewport(int x0, int y0, int w, int h) {
    /*
     * construct Viewport transformation matrix
     * assuming depth normalized as [0, 255] for visualization
     */
    Viewport = Matrix::identity(4);
    Viewport[0][3] = x0 + w / 2.f;
    Viewport[1][3] = y0 + h / 2.f;
    Viewport[2][3] = 255.f / 2.f;

    Viewport[0][0] = w / 2.f;
    Viewport[1][1] = h / 2.f;
    Viewport[2][2] = 255.f / 2.f;
}

Vec3f barycentric(Vec3f *pts, Vec2f P) {
    /*
     * calculate the barycentric coordinates given a triangle ABC and point P on a plane
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

// rasterizer, drawing triangles iterating bbox
void triangle(Vec4f *h_pts, IShader &shader, TGAImage &image, TGAImage &zbuffer) {
    /*
     * rasterize i-th face from mesh
     * h_pts assumed to be processed by vertex shader already
     */
    // normalize homogenous coords
    Vec3f pts[3];
    for (int i = 0; i < 3; ++i) {
        pts[i] = to_cartesian(h_pts[i]);
    }
    // find bounding box of a triangle in the screen plane
    Vec2f bboxmin(image.width() - 1, image.height() - 1);
    Vec2f bboxmax(0, 0);
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 2; ++j) {
            bboxmin[j] = std::max(0.f, std::min(bboxmin[j], pts[i][j]));
            bboxmax[j] = std::min(image.width() - 1.f, std::max(bboxmax[j], pts[i][j]));
        }
    }

    Vec2i P;
    // for vertex shader to fill in
    TGAColor color;
    for (P.y = bboxmin.y; P.y <= bboxmax.y; P.y++) {
        for (P.x = bboxmin.x; P.x <= bboxmax.x ; P.x++) {
            Vec3f bary_coords = barycentric(pts, P);
            // depth interpolation
            float z = bary_coords.x * pts[0].z + bary_coords.y * pts[1].z + bary_coords.z * pts[2].z;
            // truncate to uint8
            int depth = std::max(0, std::min(255, int(z + 0.5)));
            // pytorch3d frame like
            if (bary_coords.x < 0. || bary_coords.y < 0. || bary_coords.z < 0. || zbuffer.get(P.x, P.y)[0] > depth) continue;
            bool discard = shader.fragment(bary_coords, color);
            if (!discard) {
                zbuffer.set(P.x, P.y, TGAColor(depth));
                image.set(P.x, P.y, color);
            }
        }
    }
}