//
// Created by csh on 2/26/25.
// tinyrasterizer
//

#include <vector>
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

    Matrix uniform_M_T_inv;

    Vec4f vertex(int iface, int jvert) override {
        Vec4f gl_vertex = to_homogeneous(model->vert(iface, jvert));
        // vertex transformation
        gl_vertex = Viewport * Projection * ModelView * gl_vertex;
        Vec3f normal = model->normal(iface, jvert);
        normal = uniform_M_T_inv * to_homogeneous(normal);
        normal.normalize();
        // ensure cos_theta gt 0
        varying_intensity[jvert] = std::max(0.f, std::min(normal * direct_light, 1.f));
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

    Matrix uniform_M_T_inv;

    Vec4f vertex(int iface, int jvert) override {
        Vec4f gl_vertex = to_homogeneous(model->vert(iface, jvert));
        // vertex transformation
        gl_vertex = Viewport * Projection * ModelView * gl_vertex;
        // ensure cos_theta gt 0
//        Vec3f normal = model->sample_normal(iface, jvert);
        Vec3f normal = model->normal(iface, jvert);
        normal = uniform_M_T_inv * to_homogeneous(normal);
        normal.normalize();
        varying_intensity[jvert] = std::max(0.f, std::min(normal * direct_light, 1.f));
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
    Vec2f varying_uvs[3];
    Vec3f varying_normals[3];
    Vec3f varying_triangle_coords[3];

    Matrix uniform_M_T_inv;

    Vec4f vertex(int iface, int jvert) override {
        Vec4f gl_vertex = to_homogeneous(model->vert(iface, jvert));
        // vertex transformation
        gl_vertex = Viewport * Projection * ModelView * gl_vertex;
        // ensure cos_theta gt 0
        varying_uvs[jvert] = model->uv(iface, jvert);
        // to camera frame
        varying_normals[jvert] = ModelView * to_homogeneous(model->normal(iface, jvert));
        varying_triangle_coords[jvert] = ModelView * to_homogeneous(model->vert(iface, jvert));
        return gl_vertex;
    }

    bool fragment(Vec3f bary_coords, TGAColor &color) override {
        // normal interpolation
        Vec2f uv = varying_uvs[0] * bary_coords[0] + varying_uvs[1] * bary_coords[1] + varying_uvs[2] * bary_coords[2];

        Vec3f normal;
        // normal computation
        // prefer tangent space normal map
        if (model->tangent_normal().is_valid()) {
            Vec3f bary_normal = varying_normals[0] * bary_coords[0] + varying_normals[1] * bary_coords[1] + varying_normals[2] * bary_coords[2];
            Matrix A(3, 3);
            A.set_row(0, varying_triangle_coords[1] - varying_triangle_coords[0]);
            A.set_row(1, varying_triangle_coords[2] - varying_triangle_coords[0]);
            A.set_row(2, bary_normal);
            Matrix A_inv = A.inverse();
            Vec3f b_u = Vec3f(varying_uvs[1][0] - varying_uvs[0][0], varying_uvs[2][0] - varying_uvs[0][0], 0);
            Vec3f b_v = Vec3f(varying_uvs[1][1] - varying_uvs[0][1], varying_uvs[2][1] - varying_uvs[0][1], 0);
            Vec3f tangent = A_inv * b_u;
            Vec3f bitangent = A_inv * b_v;
            Matrix B(3, 3);
            B.set_col(0, tangent.normalize());
            B.set_col(1, bitangent.normalize());
            B.set_col(2, bary_normal.normalize());
            Vec3f coords = model->sample_tangent_normal(uv);
            // within camera frame
            normal = B * coords;
            normal.normalize();
        } else if (model->normal().is_valid()) {
            // then global space normal map
            normal = model->sample_normal(uv);
            // transform to camera frame
            normal = uniform_M_T_inv * to_homogeneous(normal);
            normal.normalize();
        } else {
            // within camera frame
            // use interpolated normal directly
            normal = varying_normals[0] * bary_coords[0] + varying_normals[1] * bary_coords[1] + varying_normals[2] * bary_coords[2];
            normal.normalize();
        }

        float spec_ratio = 0.f;
        float spec = 0.f;
        if (model->spec().is_valid()) {
            // heuristic
            spec_ratio = 0.4;
            // reflection direction
            Vec3f r = (2 * (normal * direct_light) * normal - direct_light).normalize();
            // v:[0, 0, 1] i.e. z axis direction of camera frame
            // if <r, v> < 0 then spec component is zero
            spec = std::pow(std::max(r.z, 0.f), model->sample_spec(uv));
        }

        float diff = std::max(0.f, std::min(normal * direct_light, 1.f));
        // final color = ambient + diffuse + component
        color = TGAColor(5, 5, 5) + model->sample_color(uv) * ((1 - spec_ratio) * diff + spec_ratio * spec);
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

    Matrix uniform_M = ModelView;
    direct_light = uniform_M * to_homogeneous(direct_light);
    direct_light.normalize();
    shader.uniform_M_T_inv = ModelView.transpose().inverse();

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
