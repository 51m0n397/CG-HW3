#ifndef _YOCTO_HAIR_H_
#define _YOCTO_HAIR_H_

// -----------------------------------------------------------------------------
// INCLUDES
// -----------------------------------------------------------------------------


#include <yocto/yocto_sampling.h>
#include <yocto/yocto_color.h>

// -----------------------------------------------------------------------------
// CONSTANTS
// -----------------------------------------------------------------------------
namespace yocto {

   static const int pMax = 3;
   static const float SqrtPiOver8 = 0.626657069f;

} // namespace yocto

// -----------------------------------------------------------------------------
// SHADING FUNCTIONS
// -----------------------------------------------------------------------------
namespace yocto {

    vec3f eval_hair(const vec3f& outgoing, const vec3f& incoming, float eta, 
                    float h, vec3f sigma_a, const vec3f& cos2kAlpha, 
                    const vec3f& sin2kAlpha, float s, const array<float, pMax + 1> v, 
                    float gammaO, const frame3f& frame);

    vec3f sample_hair(const vec3f& outgoing, const vec2f& rn, float eta, 
                      float h, vec3f sigma_a, const vec3f& cos2kAlpha, 
                      const vec3f& sin2kAlpha, float s, const array<float, pMax + 1> v, 
                      float gammaO, const frame3f& frame);

    float sample_hair_pdf(const vec3f& outgoing, const vec3f& incoming, float eta, 
                          float h, vec3f sigma_a, const vec3f& cos2kAlpha, 
                          const vec3f& sin2kAlpha, float s, const array<float, pMax + 1> v, 
                          float gammaO, const frame3f& frame);

} // namespace yocto

// -----------------------------------------------------------------------------
// SIGMA_A CONVERSION FUNCTIONS
// -----------------------------------------------------------------------------
namespace yocto {

    vec3f SigmaAFromConcentration(float ce, float cp);

    vec3f SigmaAFromReflectance(vec3f c, float beta_n);

} // namespace yocto

// -----------------------------------------------------------------------------
// UTILITY FUNCTIONS
// -----------------------------------------------------------------------------
namespace yocto {

    inline float Sqr(float v) { return v * v; }

    template <int n>
    static float Pow(float v) {
        float n2 = Pow<n / 2>(v);
        return n2 * n2 * Pow<n & 1>(v);
    }

    template <>
    inline float Pow<1>(float v) {
        return v;
    }
    template <>
    inline float Pow<0>(float v) {
        return 1;
    }

    inline float SafeASin(float x) {
        return asin(clamp(x, -1.0f, 1.0f));
    }

    inline float SafeSqrt(float x) {
        return sqrt(max(0.0f, x));
    }

    static uint32_t Compact1By1(uint32_t x) {
        // x = -f-e -d-c -b-a -9-8 -7-6 -5-4 -3-2 -1-0
        x &= 0x55555555;
        // x = --fe --dc --ba --98 --76 --54 --32 --10
        x = (x ^ (x >> 1)) & 0x33333333;
        // x = ---- fedc ---- ba98 ---- 7654 ---- 3210
        x = (x ^ (x >> 2)) & 0x0f0f0f0f;
        // x = ---- ---- fedc ba98 ---- ---- 7654 3210
        x = (x ^ (x >> 4)) & 0x00ff00ff;
        // x = ---- ---- ---- ---- fedc ba98 7654 3210
        x = (x ^ (x >> 8)) & 0x0000ffff;
        return x;
    }

    static vec2f DemuxFloat(float f) {
        uint64_t v = f * (1ull << 32);
        uint32_t bits[2] = {Compact1By1(v), Compact1By1(v >> 1)};
        return {bits[0] / float(1 << 16), bits[1] / float(1 << 16)};
    }

    inline float FrDielectric(float cosThetaI, float etaI, float etaT) {
        cosThetaI = clamp(cosThetaI, -1.0f, 1.0f);
        // Potentially swap indices of refraction
        bool entering = cosThetaI > 0.f;
        if (!entering) {
            swap(etaI, etaT);
            cosThetaI = abs(cosThetaI);
        }

        // Compute _cosThetaT_ using Snell's law
        float sinThetaI = sqrt(max((float)0, 1 - cosThetaI * cosThetaI));
        float sinThetaT = etaI / etaT * sinThetaI;

        // Handle total internal reflection
        if (sinThetaT >= 1) return 1;
        float cosThetaT = sqrt(max((float)0, 1 - sinThetaT * sinThetaT));
        float Rparl = ((etaT * cosThetaI) - (etaI * cosThetaT)) /
                    ((etaT * cosThetaI) + (etaI * cosThetaT));
        float Rperp = ((etaI * cosThetaI) - (etaT * cosThetaT)) /
                    ((etaI * cosThetaI) + (etaT * cosThetaT));
        return (Rparl * Rparl + Rperp * Rperp) / 2;
    }

    inline float AbsCosTheta(const vec3f &w) { return abs(w.z); }

} // namespace yocto

#endif