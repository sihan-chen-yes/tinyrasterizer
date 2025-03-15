//
// Created by csh on 2/26/25.
// obj file handling
//

#ifndef TINYRENDERER_MODEL_H
#define TINYRENDERER_MODEL_H

#include <vector>
#include "geometry.h"
#include "tgaimage.h"
class Model {
private:
    std::vector<Vec3f> verts_;  //xyz
    std::vector<Vec2f> uvs_;    //uv
    std::vector<Vec3f> normals_;//normal
    std::vector<std::vector<int> > faces_vert;  // per face vert idx
    std::vector<std::vector<int> > faces_uv;    // per face uv idx
    std::vector<std::vector<int> > faces_normal;// per face normal idx
    TGAImage diffusemap{};         // diffuse color texture
    TGAImage normalmap{};          // normal map
    TGAImage specmap{};         // specular map
    void load_texture(const std::string filename, const std::string suffix, TGAImage &img);
public:
    Model(const char *filename);
    ~Model();
    int nverts();
    int nfaces();
    Vec3f vert(int i);
    Vec3f vert(int iface, int nvert);
    Vec3f normal(int iface, int nvert);
    Vec2f uv(int iface, int nvert);
    std::vector<int> face(int idx);
    TGAColor sample_color(Vec2f uv);
    TGAColor sample_color(int iface, int nvert);
    Vec3f sample_normal(Vec2f uv);
    Vec3f sample_normal(int iface, int nvert);
    const TGAImage& diffuse()  const { return diffusemap;  }
    const TGAImage& normal()  const { return normalmap;  }
};


#endif //TINYRENDERER_MODEL_H
