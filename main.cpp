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

std::vector<Model*> models;
// for shadow mapping
float *shadowbuffer = NULL;
// directional light, assuming away from evaluation point
Vec3f direct_light = Vec3f(1, 1, 1).normalize();
const int width  = 800;
const int height = 800;
TGAImage ambient_occlusion(width, height, TGAImage::RGB);

// camera location
//Vec3f eye(1, 1, 3);
//Vec3f center(0, 0, 0);
//Vec3f up(0, 1, 0);

Vec3f eye(1.2,-.8,3);
Vec3f center(0,0,0);
Vec3f up(0,1,0);

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
    Vec3f varying_screen_coords[3];

    Matrix uniform_M_T_inv;
    Matrix uniform_M_shadow;

    Vec4f vertex(int iface, int jvert) override {
        Vec4f gl_vertex = to_homogeneous(model->vert(iface, jvert));
        // vertex transformation
        gl_vertex = Viewport * Projection * ModelView * gl_vertex;
        // ensure cos_theta gt 0
        varying_uvs[jvert] = model->uv(iface, jvert);
        // to camera frame
        varying_normals[jvert] = ModelView * to_homogeneous(model->normal(iface, jvert));
        varying_triangle_coords[jvert] = ModelView * to_homogeneous(model->vert(iface, jvert));
        varying_screen_coords[jvert] = Viewport * Projection * ModelView * to_homogeneous(model->vert(iface, jvert));
        return gl_vertex;
    }

    bool fragment(Vec3f bary_coords, TGAColor &color) override {
        // normal interpolation
        Vec2f uv = varying_uvs[0] * bary_coords[0] + varying_uvs[1] * bary_coords[1] + varying_uvs[2] * bary_coords[2];
        Vec3f p_screen = varying_screen_coords[0] * bary_coords[0] + varying_screen_coords[1] * bary_coords[1] + varying_screen_coords[2] * bary_coords[2];

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

        float shadow = 1;
        if (shadowbuffer != NULL) {
            //shadow mapping
            Vec3f p = varying_triangle_coords[0] * bary_coords[0] + varying_triangle_coords[1] * bary_coords[1] + varying_triangle_coords[2] * bary_coords[2];
            Vec3f p_shadow = uniform_M_shadow * to_homogeneous(p);
            int shadow_map_idx = int(p_shadow[0] + 0.5) + width * int(p_shadow[1] + 0.5);
            shadow = 0.3 + 0.7 * (shadowbuffer[shadow_map_idx] < p_shadow[2] + 50);
        }

        float spec_ratio = 0.f;
        float spec = 0.f;
        if (model->spec().is_valid()) {
            // heuristic
            spec_ratio = 0.4;
            // reflection direction
            Vec3f r = (2 * (normal * direct_light) * normal - direct_light).normalize();
            // v:[0, 0, 1] i.e. z axis direction of camera frame
            // if <r, v> lt 0 then spec component is zero
            // specular exponent should ge 1
            spec = std::pow(std::max(r.z, 0.f), std::max<float> (1, model->sample_spec(uv)));
        }
        float diff_ratio = 0.6;
        float diff = std::max(0.f, std::min(normal * direct_light, 1.f));
        // final color = ambient + diffuse + component
        // default ambient light
        TGAColor ambient_light = TGAColor (100, 100, 100);
        float amb_ratio = 0.01;
        if (ambient_occlusion.is_valid()) {
            ambient_light = ambient_occlusion.get(int(p_screen[0] + 0.5), int(p_screen[1] + 0.5));
        }
        color = ambient_light * amb_ratio + model->sample_color(uv) * (diff_ratio * diff + spec_ratio * spec) * shadow;
        // not discard
        return false;
    }
};

struct DepthShader : public IShader {
    Vec3f varying_screen_coords[3];

    Vec4f vertex(int iface, int jvert) override {
        Vec4f gl_vertex = to_homogeneous(model->vert(iface, jvert));
        // vertex transformation
        gl_vertex = Viewport * Projection * ModelView * gl_vertex;
        // ensure cos_theta gt 0
        varying_screen_coords[jvert] = Vec3f (gl_vertex);
        return gl_vertex;
    }

    virtual bool fragment(Vec3f bary_coords, TGAColor &color) {
        Vec3f p = varying_screen_coords[0] * bary_coords[0] + varying_screen_coords[1] * bary_coords[1] + varying_screen_coords[2] * bary_coords[2];
        float depth = std::max(p.z, 0.f);
        color = TGAColor(255, 255, 255) * (depth / max_depth);
        return false;
    }
};

float max_elevation_angle(float *zbuffer, Vec2f p, Vec2f dir) {
    /*
     * from p to search the height change along dir, in 2D plane
     */
    float max_angle = 0;
    for (float t = 0.; t < 1000.; t+=1.) {
        Vec2f cur = p + dir * t;
        if (cur.x >= width || cur.x < 0 || cur.y >= height || cur.y < 0) return max_angle;
        float distance = (p - cur).norm();
        // ignore starting stage
        if (distance < 1.f) continue;
        float elevation = (zbuffer[int(cur.x + 0.5) + int(cur.y + 0.5) * width] - zbuffer[int(p.x + 0.5) + int(p.y + 0.5) * width]) / max_depth;
        // slope
        max_angle = std::max(max_angle, atanf(elevation / distance));
    }
    return max_angle;
}

int main(int argc, char** argv) {
    if (2 > argc) {
        std::cerr << "Usage: " << argv[0] << " obj/model.obj" << std::endl;
        return 1;
    }

    // save all models
    for (int m = 1; m < argc; m++) {
        models.push_back(new Model(argv[m]));
    }

    // three-pass
    // first pass for ambient occlusion
    // two-pass shadow mapping
    // depth buffer initialization
    float *zbuffer = new float[width * height];
    shadowbuffer = new float[width * height];
    float *ambient_zbuffer = new float[width * height];
    // the bigger the closer to camera
    for (int i = 0; i < width * height; ++i) {
        zbuffer[i] = -std::numeric_limits<float>::max();
        shadowbuffer[i] = -std::numeric_limits<float>::max();
        ambient_zbuffer[i] = -std::numeric_limits<float>::max();
    }

    // screen based ambient occlusion pass
    // MVPV setting
    lookat(eye, center, up);
    projection(-1.f / (eye - center).norm());
    viewport(width / 8, height / 8, width * 3 / 4, height * 3 / 4);
    DepthShader ambientshader;

    rasterize(models, ambientshader, ambient_occlusion, ambient_zbuffer);

    // prepare screen space ambient occlusion
    for (int x = 0; x < width; ++x) {
        for (int y = 0; y < height; ++y) {
            // ignore background
            if (ambient_zbuffer[x + y * width] < 0) continue;
            float total = 0;
            // 8 rays, step : PI / 4, wo repeating
            for (float a = 0; a < M_PI * 2 - EPSILON; a += M_PI / 4) {
                // the bigger angle, the more hidden, the less ambient light
                total += M_PI / 2 - max_elevation_angle(ambient_zbuffer, Vec2f(x, y), Vec2f(cos(a), sin(a)));
            }
            // normalization to [0, 1]
            total /= (M_PI / 2) * 8;
            total = pow(total, 500);
            ambient_occlusion.set(x, y, TGAColor(255, 255, 255) * total);
        }
    }

    ambient_occlusion.write_tga_file("ambient_occlusion.tga");

    // first pass to generate depth map
    TGAImage depth(width, height, TGAImage::RGB);
    // MVPV setting from light source
    lookat(direct_light, center, up);
    viewport(width / 8, height / 8, width * 3 / 4, height * 3 / 4);
    projection(0);
    DepthShader depthshader;

    rasterize(models, depthshader, depth, shadowbuffer);

    depth.write_tga_file("depth.tga");
    // transformation to shadow mapping camera screen coords
    Matrix M = Viewport * Projection * ModelView;

    // second pass for rasterization
    // MVPV setting
    lookat(eye, center, up);
    projection(-1.f / (eye - center).norm());
    viewport(width / 8, height / 8, width * 3 / 4, height * 3 / 4);

    TGAImage image (width, height, TGAImage::RGB);
    // zbuffer for min depth recording
    BlinnPhongShader shader;
    shader.shadowbuffer = shadowbuffer;

    Matrix uniform_M = ModelView;
    direct_light = direct_light.normalize();
    direct_light = uniform_M * to_homogeneous(direct_light);
    direct_light.normalize();
    shader.uniform_M_T_inv = ModelView.transpose().inverse();
    shader.uniform_M_shadow = M * ModelView.inverse();

    rasterize(models, shader, image, zbuffer);
    image.write_tga_file("shader.tga");

    for (Model *model: models) {
        delete model;
    }
    delete [] shadowbuffer;
    delete [] zbuffer;
    delete [] ambient_zbuffer;
    return 0;
}
