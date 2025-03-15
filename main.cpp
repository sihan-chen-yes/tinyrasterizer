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
#include "our_gl.h"

Model *model = NULL;
// directional light, assuming away from evaluation point
Vec3f direct_light = Vec3f(1, 1, 1).normalize();
const int width  = 800;
const int height = 800;
// camera location
Vec3f eye(1, 1, 3);
Vec3f center(0, 0, 0);
Vec3f up(0, 1, 0);

struct ToonShader: public IShader {
    /*
     * refresh for each triangle face
     * written by vertex shader, read by fragment shader
     */
    Vec3f varying_intensity;

    Vec4f vertex(int iface, int jvert) override {
        Vec4f gl_vertex = to_homogeneous(model->vert(iface, jvert));
        // vertex transformation
        gl_vertex = Viewport * Projection * ModelView * gl_vertex;
        // ensure cos_theta gt 0
        varying_intensity[jvert] = std::max(0.f, model->normal(iface, jvert) * direct_light);
        return gl_vertex;
    }

    bool fragment(Vec3f bary_coords, TGAColor &color) override {
        // Gouraud intensity interpolation
        float intensity = varying_intensity * bary_coords;
        // fixed intensity levels
        if (intensity > 0.85) intensity = 1;
        else if (intensity > 0.60) intensity = 0.80;
        else if (intensity > 0.45) intensity = 0.60;
        else if (intensity > 0.30) intensity = 0.45;
        else if (intensity > 0.15) intensity = 0.30;
        else intensity = 0;
        color = TGAColor(255, 155, 0) * intensity;
        // not discard
        return false;
    }
};


struct GouraudShader: public IShader {
    /*
     * refresh for each triangle face
     * written by vertex shader, read by fragment shader
     */
    Vec3f varying_intensity;
    Vec2f varying_uvs[3];

    Vec4f vertex(int iface, int jvert) override {
        Vec4f gl_vertex = to_homogeneous(model->vert(iface, jvert));
        // vertex transformation
        gl_vertex = Viewport * Projection * ModelView * gl_vertex;
        // ensure cos_theta gt 0
        varying_intensity[jvert] = std::max(0.f, model->normal(iface, jvert) * direct_light);
        varying_uvs[jvert] = model->uv(iface, jvert);
        return gl_vertex;
    }

    bool fragment(Vec3f bary_coords, TGAColor &color) override {
        // Gouraud intensity interpolation
        float intensity = varying_intensity * bary_coords;
        Vec2f uv = varying_uvs[0] * bary_coords[0] + varying_uvs[1] * bary_coords[1] + varying_uvs[2] * bary_coords[2];
        color = model->sample_color(uv) * intensity;
        // not discard
        return false;
    }
};

struct BlinnPhongShader: public IShader {
    /*
     * refresh for each triangle face
     * written by vertex shader, read by fragment shader
     */
    Vec3f varying_normals[3];
    Vec2f varying_uvs[3];

    Vec4f vertex(int iface, int jvert) override {
        Vec4f gl_vertex = to_homogeneous(model->vert(iface, jvert));
        // vertex transformation
        gl_vertex = Viewport * Projection * ModelView * gl_vertex;
        // ensure cos_theta gt 0
        varying_normals[jvert] = model->normal(iface, jvert);
        varying_uvs[jvert] = model->uv(iface, jvert);
        return gl_vertex;
    }

    bool fragment(Vec3f bary_coords, TGAColor &color) override {
        // normal interpolation
        Vec3f normal = varying_normals[0] * bary_coords[0] + varying_normals[1] * bary_coords[1] + varying_normals[2] * bary_coords[2];
        Vec2f uv = varying_uvs[0] * bary_coords[0] + varying_uvs[1] * bary_coords[1] + varying_uvs[2] * bary_coords[2];
        float intensity = std::max(0.f, normal * direct_light);
        color = model->sample_color(uv) * intensity;
        // not discard
        return false;
    }
};

int main(int argc, char** argv) {
    if (2 == argc) {
        model = new Model(argv[1]);
    } else {
        model = new Model("obj/african_head.obj");
    }

    // MVPV setting
    lookat(eye, center, up);
    projection(-1.f / (eye - center).norm());
    viewport(width / 8, height / 8, width * 3 / 4, height * 3 / 4);

    TGAImage image  (width, height, TGAImage::RGB);
    // zbuffer for min depth recording
    TGAImage zbuffer(width, height, TGAImage::GRAYSCALE);
    BlinnPhongShader shader;

    for (int i = 0; i < model->nfaces(); ++i) {
        Vec4f screen_coords[3];
        for (int j = 0; j < 3; ++j) {
            screen_coords[j] = shader.vertex(i, j);
        }
        triangle(screen_coords, shader, image, zbuffer);
    }

    image.write_tga_file("shader_test.tga");
    zbuffer.write_tga_file("zbuffer_debug.tga");

    delete model;
    return 0;
}
