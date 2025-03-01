//
// Created by csh on 2/26/25.
// tga image handling
//
#ifndef TINYRENDERER_TGAIMAGE_H
#define TINYRENDERER_TGAIMAGE_H
#pragma once
#include <cstdint>
#include <fstream>
#include <vector>
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
    TGAColor operator *(float f)  { return TGAColor(bgra[2] * f, bgra[1] * f, bgra[0] * f); }
    TGAColor operator +(TGAColor c)  { return TGAColor(bgra[2] + c[2], bgra[1] + c[1], bgra[0] + c[0]); }

    TGAColor() {}
    TGAColor(std::uint8_t r, std::uint8_t g, std::uint8_t b, std::uint8_t a = 255, std::uint8_t bpp = 4)
            : bgra{b, g, r, a}, bytespp(bpp) {
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
private:
    bool   load_rle_data(std::ifstream &in); //decompression
    bool unload_rle_data(std::ofstream &out) const; //compression
    int w = 0, h = 0;
    std::uint8_t bpp = 0;
    std::vector<std::uint8_t> data = {};
};
#endif //TINYRENDERER_TGAIMAGE_H
