//
// Created by csh on 2/26/25.
// obj file handling
//

#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <vector>
#include "model.h"

Model::Model(const char *filename) : verts_(), faces_vert() {
    std::ifstream in;
    in.open (filename, std::ifstream::in);
    if (in.fail()) return;
    std::string line;
    while (!in.eof()) {
        std::getline(in, line);
        std::istringstream iss(line.c_str());
        char trash;
        // if equal return 0
        if (!line.compare(0, 2, "v ")) {
            iss >> trash;
            Vec3f v;
            for (int i=0;i<3;i++) iss >> v.raw[i];
            verts_.push_back(v);
        } else if (!line.compare(0, 3, "vn ")) {
            iss >> trash >> trash;
            Vec3f n;
            for (int i=0;i<3;i++) iss >> n.raw[i];
            n.normalize();
            normals_.push_back(n);
        } else if (!line.compare(0, 3, "vt ")) {
            iss >> trash >> trash;
            Vec2f uv;
            for (int i=0;i<2;i++) iss >> uv.raw[i];
            uvs_.push_back({uv.x, uv.y});
        } else if (!line.compare(0, 2, "f ")) {
            std::vector<int> face_vert, face_uv, face_normal;
            int itrash, f, t, n;
            iss >> trash;
            int cnt = 0;
            while (iss >> f >> trash >> t >> trash >> n) {
                // in wavefront obj all indices start at 1, not zero
                face_vert.push_back(--f);
                face_uv.push_back(--t);
                face_normal.push_back(--n);
                cnt++;
            }
            if (3!=cnt) {
                std::cerr << "Error: the obj file is supposed to be triangulated" << std::endl;
                return;
            }
            faces_vert.push_back(face_vert);
            faces_uv.push_back(face_uv);
            faces_normal.push_back(face_normal);
        }
    }
    std::cerr << "# v# " << nverts() << " f# "  << nfaces() << " vt# " << uvs_.size() << " vn# " << normals_.size() << std::endl;
    load_texture(filename, "_diffuse.tga", diffusemap);
    load_texture(filename, "_nm.tga", normalmap);
    load_texture(filename, "_spec.tga", specmap);
}

Model::~Model() {
}

int Model::nverts() {
    return (int)verts_.size();
}

int Model::nfaces() {
    return (int)faces_vert.size();
}

std::vector<int> Model::face(int idx) {
    return faces_vert[idx];
}

Vec3f Model::vert(int i) {
    return verts_[i];
}

Vec3f Model::vert(int iface, int jvert) {
    return verts_[faces_vert[iface][jvert]];
}

Vec3f Model::normal(int iface, int jvert) {
    return normals_[faces_normal[iface][jvert]];
}

Vec2f Model::uv(int iface, int jvert) {
    return uvs_[faces_uv[iface][jvert]];
}

TGAColor Model::sample_color(Vec2f uv) {
    return diffusemap.get(uv.u * diffusemap.width(), uv.v * diffusemap.height());
}

TGAColor Model::sample_color(int iface, int jvert) {
    return sample_color(uv(iface, jvert));
}

Vec3f Model::sample_normal(Vec2f uv) {
    TGAColor c = normalmap.get(uv.u * diffusemap.width(), uv.v * diffusemap.height());
    return Vec3f (c[2], c[1], c[0]) * 2.f / 255.f - Vec3f (1, 1, 1);
}

Vec3f Model::sample_normal(int iface, int jvert) {
    return sample_normal(uv(iface, jvert));
}

void Model::load_texture(std::string filename, const std::string suffix, TGAImage &img) {
    size_t dot = filename.find_last_of(".");
    if (dot==std::string::npos) return;
    std::string texfile = filename.substr(0,dot) + suffix;
    std::cerr << "texture file " << texfile << " loading " << (img.read_tga_file(texfile.c_str()) ? "ok" : "failed") << std::endl;
    img.flip_vertically();
}