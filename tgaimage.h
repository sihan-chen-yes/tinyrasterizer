//
// Created by csh on 2/26/25.
// tga image handling
//
#ifndef TINYRENDERER_TGAIMAGE_H
#define TINYRENDERER_TGAIMAGE_H
#pragma once
#include <cstdint>
#include <fstream>
#include <iostream>
#include <vector>
#include "geometry.h"
// TGA picture rules, can't pad extra bytes
#pragma pack(push,1)
struct TGAHeader {
    std::uint8_t  idlength = 0;
    std::uint8_t  colormaptype = 0;
    std::uint8_t  datatypecode = 0;
    std::uint16_t colormaporigin = 0;
    std::uint16_t colormaplength = 0;
    std::uint8_t  colormapdepth = 0;
    std::uint16_t x_origin = 0;
    std::uint16_t y_origin = 0;
    std::uint16_t width = 0;
    std::uint16_t height = 0;
    std::uint8_t  bitsperpixel = 0;
    std::uint8_t  imagedescriptor = 0;
};
#pragma pack(pop)

struct TGAColor {
    std::uint8_t bgra[4] = {0, 0, 0, 0}; // BGRA order
    std::uint8_t bytespp = 4; // bytes per pixel
    std::uint8_t& operator[](const int i) { return bgra[i]; }

    TGAColor operator *(float f)  {
        float R = bgra[2] * f;
        float G = bgra[1] * f;
        float B = bgra[0] * f;

        if (valid_range(R, G, B)) {
            return TGAColor(R, G, B, bgra[3]);
        } else {
            std::cerr << "[Warn] color clamped to [0..255] in operator*()\n";
            R = std::clamp<float>(R, 0, 255);
            G = std::clamp<float>(G, 0, 255);
            B = std::clamp<float>(B, 0, 255);
            return TGAColor(R, G, B, bgra[3]);
        }
    }

    TGAColor operator *(Vec3f f)  {
        float R = bgra[2] * (f[0] + f[1] + f[2]);
        float G = bgra[1] * (f[0] + f[1] + f[2]);
        float B = bgra[0] * (f[0] + f[1] + f[2]);

        if (valid_range(R, G, B)) {
            return TGAColor(R, G, B, bgra[3]);
        } else {
            std::cerr << "[Warn] color clamped to [0..255] in operator*()\n";
            R = std::clamp<float>(R, 0, 255);
            G = std::clamp<float>(G, 0, 255);
            B = std::clamp<float>(B, 0, 255);
            return TGAColor(R, G, B, bgra[3]);
        }
    }

    TGAColor operator +(TGAColor c)  {
        int R = bgra[2] + c.bgra[2];
        int G = bgra[1] + c.bgra[1];
        int B = bgra[0] + c.bgra[0];

        if (valid_range(R, G, B)) {
            return TGAColor(R, G, B, bgra[3]);
        } else {
            std::cerr << "[Warn] color clamped to 255 in operator+()\n";
            R = std::clamp(R, 0, 255);
            G = std::clamp(G, 0, 255);
            B = std::clamp(B, 0, 255);
            return TGAColor(R, G, B, bgra[3]);
        }
    }

    bool valid_range(int R, int G, int B) {
        return (R >= 0 && R <= 255)
               && (G >= 0 && G <= 255)
               && (B >= 0 && B <= 255);
    }

    bool valid_range(float R, float G, float B) {
        return (R >= 0.f && R <= 255.f)
               && (G >= 0.f && G <= 255.f)
               && (B >= 0.f && B <= 255.f);
    }

    friend std::ostream& operator<<(std::ostream& s, TGAColor &c) {
        s << "R:" << int(c[2]) << ", " << "G:" << int(c[1]) << ", " << "B:" << int(c[2]) << ", " << "alpha:" << int(c[3]) << "\n";
        return s;
    }

    TGAColor() {}
    TGAColor(std::uint8_t r, std::uint8_t g, std::uint8_t b, std::uint8_t a = 255, std::uint8_t bpp = 4)
            : bgra{b, g, r, a}, bytespp(bpp) {
    }
    TGAColor(unsigned char v) : bgra(), bytespp(1) {
        for (int i=0; i<4; i++) bgra[i] = 0;
        bgra[0] = v;
    }

};

struct TGAImage {
    enum Format { GRAYSCALE=1, RGB=3, RGBA=4 };
    TGAImage() = default;
    TGAImage(const int w, const int h, const int bpp);
    bool  read_tga_file(const std::string filename);
    bool write_tga_file(const std::string filename, const bool vflip=true, const bool rle=true) const; // whether use run-length encoding for compression
    void flip_horizontally();
    void flip_vertically();
    TGAColor get(const int x, const int y) const;
    void set(const int x, const int y, const TGAColor &c);
    int width()  const;
    int height() const;
    bool is_valid() const;
private:
    bool   load_rle_data(std::ifstream &in); //decompression
    bool unload_rle_data(std::ofstream &out) const; //compression
    int w = 0, h = 0;
    std::uint8_t bpp = 0;
    std::vector<std::uint8_t> data = {};
};
#endif //TINYRENDERER_TGAIMAGE_H
