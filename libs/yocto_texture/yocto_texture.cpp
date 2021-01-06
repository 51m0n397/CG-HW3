#include "yocto_texture.h"

// -----------------------------------------------------------------------------
// LOCAL FUNCTIONS
// -----------------------------------------------------------------------------

namespace yocto {

    float erf(float x) {
        float a1 =  0.254829592;
        float a2 = -0.284496736;
        float a3 =  1.421413741;
        float a4 = -1.453152027;
        float a5 =  1.061405429;
        float p  =  0.3275911;
        
        float sign = x < 0 ? -1 : 1;
        x = abs(x);
        
        float t = 1.0/(1.0 + p*x);
        float y = 1.0 - ((((a5*t + a4)*t + a3)*t + a2)*t + a1)*t*exp(-x*x);
        
        return sign*y;
    }

    float derf(float x) {
        return 2.0 / sqrt(pi) * exp(-x * x);
    }

    float erfInv(float x) {
        float y = 0.0, err;
        do {
            err = erf(y) - x;
            y -= err/derf(y);
        } while (abs(err) > 1e-6);
        return y;
    }

    float C(float sigma) {
        return 1.0 / erf(0.5 / (sigma * sqrt(2.0)));
    }

    float truncCdfInv(float x, float sigma) {
        return 0.5 + sqrt(2) * sigma * erfInv((2.0 * x - 1.0) / C(sigma));
    }

    vec3f YCbCrToRGB(vec3f YCbCr) {
        mat3f conv = mat3f{
            {298,  298, 298},
            {  0, -100, 516},
            {408, -208,   0}};
        
        return conv*YCbCr/256.0 + vec3f{-222, 135, -276}/256.0;
    }

    float softClipContrast(float x, float W) {
        float u = x > 0.5 ? 1.0 - x : x;
        float result;
        if (u >= 0.5 - 0.25*W)
            result = (u - 0.5)/W + 0.5;
        else if (W >= 2.0/3.0)
            result = 8.0*(1.0/W - 1.0)*(u/(2.0 - W))*(u/(2.0 - W)) + (3.0 - 2.0/W)*u/(2.0 - W);
        else if (u >= 0.5 - 0.75*W)
            result = ((u - (0.5 - 0.75*W))/W)*((u - (0.5 - 0.75*W))/W);
        else
            result = 0.0;

        if (x > 0.5)
            result = 1.0 - result;

        return result;
    }

    vec3f sampleTile(const image<vec4b>& pattern, vec2f vertex, vec2f offset) {
        vec2i ij = {int(vertex.x), int(vertex.y)};
        vec4f rands = randTex[ij];
        vec2f pos = 0.25 + vec2f{rands.x, rands.y}*0.5 + offset;
        vec2i pos_i = {int(pos.x*pattern.imsize().x), int(pos.y*pattern.imsize().y)};
    return byte_to_float(xyz(pattern[pos_i]));
}

} // namespace yocto

// -----------------------------------------------------------------------------
// Randomized Texture Tiling IMPLEMENTATION
// -----------------------------------------------------------------------------

namespace yocto {

    void gaussianize_texture(pathtrace_texture* texture) {

        //converting from hdr to ldr
        if (texture->ldr.empty()) {
            texture->ldr = float_to_byte(rgb_to_srgb(texture->hdr));
            texture->hdr = {};
        }

        texture->randomize = true;

        auto minSize = min(texture->ldr.imsize());
        auto size = 1;
        while (2*size <= minSize)
            size *= 2;
        auto w = size;
        auto h = size;

        // Converting to ycbcr and computing histogram
        int histo[256];
        for (auto i = 0; i < 256; ++i)
        histo[i] = 0;

        for (auto i = 0; i < w*h; ++i) {
            auto o = rgbToYcbcrOff;
            auto M = rgbToYcbcrMat;
            auto r = texture->ldr[i].x;
            auto g = texture->ldr[i].y;
            auto b = texture->ldr[i].z;

            auto Y  = o[0] + M[0][0]*r + M[0][1]*g + M[0][2]*b;
            auto Cb = o[1] + M[1][0]*r + M[1][1]*g + M[1][2]*b;
            auto Cr = o[2] + M[2][0]*r + M[2][1]*g + M[2][2]*b;
            texture->ldr[i].x = byte(max(0.0f, min(255.0f, floor(Y ))));
            texture->ldr[i].y = byte(max(0.0f, min(255.0f, floor(Cb))));
            texture->ldr[i].z = byte(max(0.0f, min(255.0f, floor(Cr))));
            
            histo[texture->ldr[i].x]++;
        }

        for (auto i = 1; i < 256; ++i)
            histo[i] += histo[i - 1];

        // building lut
        float lutF[256];
        byte lut[256];
        for (auto i = 0; i < 256; ++i) {
            lutF[i] = truncCdfInv(float(histo[i])/float(histo[255]), 1.0/6.0);
            lut[i] = max(0, min(255, (int) floor(255*lutF[i])));
        }

        // applying lut to texture
        for (auto i = 0; i < w*h; ++i)
            texture->ldr[i].x = lut[texture->ldr[i].x];

        // calculating inverse lut
        for (auto i = 0; i < ilutRes; ++i) {
            auto f = (i + 0.5)/ilutRes;
            texture->inv_lut[i] = 255;
            for (auto j = 0; j < 256; ++j) {
                if (f < lutF[j]) {
                    texture->inv_lut[i] = j;
                    break;
                }
            }
        }

    }

    vec4f lookup_randomized_texture(const pathtrace_texture* texture, const vec2f& coord) {
        vec2f pos = (coord - vec2f{0.5, 0.5});
        vec2f lattice = worldToLattice * pos;
        vec2f cell = floor(lattice);
        vec2f uv = lattice - cell;

        vec2f v0 = cell;
        if (uv.x + uv.y >= 1.0) {
            v0 += 1.0;
            uv = 1.0 - vec2f{uv.y, uv.x};
        }
        vec2f v1 = cell + vec2f{1, 0};
        vec2f v2 = cell + vec2f{0, 1};

        vec3f color0 = sampleTile(texture->ldr, v0, pos - latticeToWorld*v0);
        vec3f color1 = sampleTile(texture->ldr, v1, pos - latticeToWorld*v1);
        vec3f color2 = sampleTile(texture->ldr, v2, pos - latticeToWorld*v2);

        vec3f uvw = vec3f{1.0f - uv.x - uv.y, uv.x, uv.y};
        uvw = uvw * uvw * uvw;
        uvw /= uvw.x + uvw.y + uvw.z;

        vec3f YCbCr = uvw.x*color0 + uvw.y*color1 + uvw.z*color2;

        YCbCr.x = softClipContrast(YCbCr.x, uvw.x + uvw.y + uvw.z);
        YCbCr.x = byte_to_float(texture->inv_lut[int(YCbCr.x*1024)]);
        
        vec3f color = srgb_to_rgb(YCbCrToRGB(YCbCr));

        return {color.x, color.y, color.z, 1.0};
    }

} // namespace yocto