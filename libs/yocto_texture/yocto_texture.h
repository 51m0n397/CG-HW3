#ifndef _YOCTO_TEXTURE_H_
#define _YOCTO_TEXTURE_H_

// -----------------------------------------------------------------------------
// INCLUDES
// -----------------------------------------------------------------------------

#include <yocto/yocto_math.h>
#include <yocto/yocto_image.h>
#include <yocto/yocto_color.h>
#include <yocto/yocto_sampling.h>
#include <yocto_pathtrace/yocto_pathtrace.h>

// -----------------------------------------------------------------------------
// UTILITY FUNCTIONS
// -----------------------------------------------------------------------------

namespace yocto {

    inline image<vec4f> make_random_tex(int res) {
        auto size = vec2i{res, res};
        auto randImg = image<vec4f>{size};
        auto rng = make_rng(1301081);
        
        for (int i=0; i<res; i++) {
            for (int j=0; j<res; j++) {
                randImg[{i, j}] = {rand1f(rng), rand1f(rng), rand1f(rng), rand1f(rng)};
            }
        }

        return randImg;
    };

    inline float floor(float a) { return std::floorf(a); }

    inline vec2f floor(const vec2f& a) { return {floor(a.x), floor(a.y)}; }

    inline vec3f floor(const vec3f& a) { return {floor(a.x), floor(a.y), floor(a.z)}; }

} // namespace yocto


// -----------------------------------------------------------------------------
// CONSTANTS
// -----------------------------------------------------------------------------
namespace yocto {

    const auto rgbToYcbcrOff = vec3f{16, 128, 128};

    const auto rgbToYcbcrMat = mat3f{
        {  65.5/255.0,  128.6/255.0,   25.0/255.0},
        {- 37.8/255.0, - 74.2/255.0,  112.0/255.0},
        { 112.0/255.0, - 93.8/255.0, - 18.2/255.0}
    };

    const auto latticeToWorld = mat2f{
        {(float) 0.25 * cos(pi / 3.0), 0.25}, 
        {(float) 0.25 * sin(pi / 3.0), 0.0 }
    };

    const auto worldToLattice = inverse(latticeToWorld);

    const auto randRes = 512;
    
    const auto randTex = make_random_tex(randRes);

    const auto ilutRes = 1024;

} // namespace yocto

// -----------------------------------------------------------------------------
// Randomized Texture Tiling
// -----------------------------------------------------------------------------

namespace yocto {

    void gaussianize_texture(pathtrace_texture* texture);

    vec4f lookup_randomized_texture(const pathtrace_texture* texture, const vec2f& coord);

} // namespace yocto

#endif