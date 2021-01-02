#include "yocto_hair.h"

// -----------------------------------------------------------------------------
// HAIR LOCAL FUNCTIONS
// -----------------------------------------------------------------------------
namespace yocto {

    float I0(float x) {
        float val = 0;
        float x2i = 1;
        int64_t ifact = 1;
        int i4 = 1;

        for (int i = 0; i < 10; ++i) {
            if (i > 1) ifact *= i;
            val += x2i / (i4 * Sqr(ifact));
            x2i *= x * x;
            i4 *= 4;
        }
        return val;
    }

    float LogI0(float x) {
        if (x > 12)
            return x + 0.5 * (-log(2 * pi) + log(1 / x) + 1 / (8 * x));
        else
            return log(I0(x));
    }

    static array<vec3f, pMax + 1> Ap(float cosThetaO, float eta, float h,
                                     const vec3f &T) {
        array<vec3f, pMax + 1> ap;
        // Compute p=0 attenuation at initial cylinder intersection
        float cosGammaO = SafeSqrt(1 - h * h);
        float cosTheta = cosThetaO * cosGammaO;
        float f = FrDielectric(cosTheta, 1.f, eta); //!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ap[0] = vec3f{f};

        // Compute p=1 attenuation term
        ap[1] = Sqr(1 - f) * T;

        // Compute attenuation terms up to p=pMax
        for (int p = 2; p < pMax; ++p) ap[p] = ap[p - 1] * T * f;

        // Compute attenuation term accounting for remaining orders of scattering
        ap[pMax] = ap[pMax - 1] * f * T / (vec3f{1.f} - T * f);
        return ap;
    }

    static float Mp(float cosThetaI, float cosThetaO, float sinThetaI,
                    float sinThetaO, float v) {
        float a = cosThetaI * cosThetaO / v;
        float b = sinThetaI * sinThetaO / v;
        float mp =
            (v <= .1)
                ? (exp(LogI0(a) - b - 1 / v + 0.6931f + log(1 / (2 * v))))
                : (exp(-b) * I0(a)) / (sinh(1 / v) * 2 * v);
        return mp;
    }

    float Phi(int p, float gammaO, float gammaT) {
        return 2 * p * gammaT - 2 * gammaO + p * pi;
    }

    float Logistic(float x, float s) {
        x = abs(x);
        return exp(-x / s) / (s * Sqr(1 + exp(-x / s)));
    }

    float LogisticCDF(float x, float s) {
        return 1 / (1 + exp(-x / s));
    }


    float TrimmedLogistic(float x, float s, float a, float b) {
        return Logistic(x, s) / (LogisticCDF(b, s) - LogisticCDF(a, s));
    }

    static float SampleTrimmedLogistic(float u, float s, float a, float b) {
        float k = LogisticCDF(b, s) - LogisticCDF(a, s);
        float x = -s * log(1 / (u * k + LogisticCDF(a, s)) - 1);
        return clamp(x, a, b);
    }

    float Np(float phi, int p, float s, float gammaO, float gammaT) {
        float dphi = phi - Phi(p, gammaO, gammaT);
        while (dphi > pi) dphi -= 2 * pi;
        while (dphi < -pi) dphi += 2 * pi;
        return TrimmedLogistic(dphi, s, -pi, pi);
    }

    array<float, pMax + 1> ComputeApPdf(float cosThetaO, float eta, float h,
                                        vec3f& sigma_a) {
        // Compute array of Ap values for cosThetaO
        float sinThetaO = SafeSqrt(1 - cosThetaO * cosThetaO);

        // Compute cosThetaT for refracted ray
        float sinThetaT = sinThetaO / eta;
        float cosThetaT = SafeSqrt(1 - Sqr(sinThetaT));

        // Compute gammaT for refracted ray
        float etap = sqrt(eta * eta - Sqr(sinThetaO)) / cosThetaO;
        float sinGammaT = h / etap;
        float cosGammaT = SafeSqrt(1 - Sqr(sinGammaT));

        // Compute the transmittance T of a single path through the cylinder
        vec3f T = exp(-sigma_a * (2 * cosGammaT / cosThetaT));
        array<vec3f, pMax + 1> ap = Ap(cosThetaO, eta, h, T);

        // Compute Ap PDF from individual Ap terms
        array<float, pMax + 1> apPdf;
        
        float sumY = 0;
        for (int i = 0; i <= pMax; ++i) sumY += luminance(ap[i]);
        for (int i = 0; i <= pMax; ++i) apPdf[i] = luminance(ap[i]) / sumY;

        return apPdf;
    }

}  // namespace yocto

// -----------------------------------------------------------------------------
// IMPLEMENTATION OF SHADING FUNCTIONS
// -----------------------------------------------------------------------------
namespace yocto {

    vec3f eval_hair(const vec3f& outgoing, const vec3f& incoming, float eta, 
                    float h, vec3f sigma_a, const vec3f& cos2kAlpha, 
                    const vec3f& sin2kAlpha, float s, const array<float, pMax + 1> v, float gammaO, const frame3f& frame) {

        vec3f wo = transform_direction(frame, outgoing);
        vec3f wi = transform_direction(frame, incoming);

        // Compute hair coordinate system terms related to wo
        float sinThetaO = wo.x;
        float cosThetaO = SafeSqrt(1 - Sqr(sinThetaO));
        float phiO = atan2(wo.z, wo.y);

        // Compute hair coordinate system terms related to wi
        float sinThetaI = wi.x;
        float cosThetaI = SafeSqrt(1 - Sqr(sinThetaI));
        float phiI = atan2(wi.z, wi.y);

         // Compute cosThetaT for refracted ray
        float sinThetaT = sinThetaO / eta;
        float cosThetaT = SafeSqrt(1 - Sqr(sinThetaT));

        // Compute gammaT for refracted ray
        float etap = sqrt(eta * eta - Sqr(sinThetaO)) / cosThetaO;
        float sinGammaT = h / etap;
        float cosGammaT = SafeSqrt(1 - Sqr(sinGammaT));
        float gammaT = SafeASin(sinGammaT);

        // Compute the transmittance _T_ of a single path through the cylinder
        vec3f T = exp(-sigma_a * (2 * cosGammaT / cosThetaT));

        // Evaluate hair BSDF
        float phi = phiI - phiO;

        array<vec3f, pMax + 1> ap = Ap(cosThetaO, eta, h, T);

        vec3f fsum = zero3f;

        for (int p = 0; p < pMax; ++p) {
            // Compute sinThetaOp and cosThetaOp terms accounting for scales
            float sinThetaOp, cosThetaOp;

            if (p == 0) {
                sinThetaOp = sinThetaO * cos2kAlpha[1] - cosThetaO * sin2kAlpha[1];
                cosThetaOp = cosThetaO * cos2kAlpha[1] + sinThetaO * sin2kAlpha[1];
            } else if (p == 1) {
                sinThetaOp = sinThetaO * cos2kAlpha[0] + cosThetaO * sin2kAlpha[0];
                cosThetaOp = cosThetaO * cos2kAlpha[0] - sinThetaO * sin2kAlpha[0];
            } else if (p == 2) {
                sinThetaOp = sinThetaO * cos2kAlpha[2] + cosThetaO * sin2kAlpha[2];
                cosThetaOp = cosThetaO * cos2kAlpha[2] - sinThetaO * sin2kAlpha[2];
            } else {
                sinThetaOp = sinThetaO;
                cosThetaOp = cosThetaO;
            }

            // Handle out-of-range cosThetaOp from scale adjustment
            cosThetaOp = abs(cosThetaOp);

            fsum += Mp(cosThetaI, cosThetaOp, sinThetaI, sinThetaOp, v[p]) * ap[p] *
                    Np(phi, p, s, gammaO, gammaT);
        }

        // Compute contribution of remaining terms after pMax
        fsum += Mp(cosThetaI, cosThetaO, sinThetaI, sinThetaO, v[pMax]) * ap[pMax] /
                (2.f * pi);
        if (AbsCosTheta(wi) > 0) fsum /= AbsCosTheta(wi); //!!!!!!!!!!!!!!!!

        return fsum;
    }

    vec3f sample_hair(const vec3f& outgoing, const vec2f& rn, float eta, 
                      float h, vec3f sigma_a, const vec3f& cos2kAlpha, 
                      const vec3f& sin2kAlpha, float s, const array<float, pMax + 1> v, float gammaO, const frame3f& frame) {
        
        vec3f wo = transform_direction(frame, outgoing);
        
        // Compute hair coordinate system terms related to wo
        float sinThetaO = wo.x;
        float cosThetaO = SafeSqrt(1 - Sqr(sinThetaO));
        float phiO = atan2(wo.z, wo.y);

        // Derive four random samples from rn
        vec2f u[2] = {DemuxFloat(rn[0]), DemuxFloat(rn[1])};

        // Determine which term p to sample for hair scattering
        array<float, pMax + 1> apPdf = ComputeApPdf(cosThetaO, eta, h, sigma_a);
        int p;
        for (p = 0; p < pMax; ++p) {
            if (u[0][0] < apPdf[p]) break;
            u[0][0] -= apPdf[p];
        }

        // Rotate sinThetaO and cosThetaO to account for hair scale tilt
        float sinThetaOp, cosThetaOp;
        if (p == 0) {
            sinThetaOp = sinThetaO * cos2kAlpha[1] - cosThetaO * sin2kAlpha[1];
            cosThetaOp = cosThetaO * cos2kAlpha[1] + sinThetaO * sin2kAlpha[1];
        }
        else if (p == 1) {
            sinThetaOp = sinThetaO * cos2kAlpha[0] + cosThetaO * sin2kAlpha[0];
            cosThetaOp = cosThetaO * cos2kAlpha[0] - sinThetaO * sin2kAlpha[0];
        } else if (p == 2) {
            sinThetaOp = sinThetaO * cos2kAlpha[2] + cosThetaO * sin2kAlpha[2];
            cosThetaOp = cosThetaO * cos2kAlpha[2] - sinThetaO * sin2kAlpha[2];
        } else {
            sinThetaOp = sinThetaO;
            cosThetaOp = cosThetaO;
        }

        // Sample Mp to compute cosThetaI
        u[1][0] = max(u[1][0], 1e-5);
        float cosTheta =
            1 + v[p] * log(u[1][0] + (1 - u[1][0]) * exp(-2 / v[p]));
        float sinTheta = SafeSqrt(1 - Sqr(cosTheta));
        float cosPhi = cos(2 * pi * u[1][1]);
        float sinThetaI = -cosTheta * sinThetaOp + sinTheta * cosPhi * cosThetaOp;
        float cosThetaI = SafeSqrt(1 - Sqr(sinThetaI));

         // Sample Np to compute dphi

        // Compute gammaT for refracted ray
        float etap = sqrt(eta * eta - Sqr(sinThetaO)) / cosThetaO;
        float sinGammaT = h / etap;
        float gammaT = SafeASin(sinGammaT);
        float dphi;
        if (p < pMax)
            dphi =
                Phi(p, gammaO, gammaT) + SampleTrimmedLogistic(u[0][1], s, -pi, pi);
        else
            dphi = 2 * pi * u[0][1];

        // Compute incoming from sampled hair scattering angles
        float phiI = phiO + dphi;
        vec3f wi = {sinThetaI, cosThetaI * cos(phiI), cosThetaI * sin(phiI)};

        vec3f incoming = transform_direction(inverse(frame), wi);

        return incoming;
    }

    float sample_hair_pdf(const vec3f& outgoing, const vec3f& incoming, float eta, 
                          float h, vec3f sigma_a, const vec3f& cos2kAlpha, 
                          const vec3f& sin2kAlpha, float s, const array<float, pMax + 1> v, float gammaO, const frame3f& frame) {
        
        vec3f wo = transform_direction(frame, outgoing);
        vec3f wi = transform_direction(frame, incoming);
        
        // Compute hair coordinate system terms related to wo
        float sinThetaO = wo.x;
        float cosThetaO = SafeSqrt(1 - Sqr(sinThetaO));
        float phiO = atan2(wo.z, wo.y);

        // Compute hair coordinate system terms related to wi
        float sinThetaI = wi.x;
        float cosThetaI = SafeSqrt(1 - Sqr(sinThetaI));
        float phiI = atan2(wi.z, wi.y);

        // Compute gammaT for refracted ray
        float etap = sqrt(eta * eta - Sqr(sinThetaO)) / cosThetaO;
        float sinGammaT = h / etap;
        float gammaT = SafeASin(sinGammaT);

        // Compute PDF for Ap terms
        array<float, pMax + 1> apPdf = ComputeApPdf(cosThetaO, eta, h, sigma_a);

        // Compute PDF sum for hair scattering events
        float phi = phiI - phiO;
        float pdf = 0;
        for (int p = 0; p < pMax; ++p) {
            // Compute sinThetaOp and cosThetaOp terms accounting for scales
            float sinThetaOp, cosThetaOp;
            if (p == 0) {
                sinThetaOp = sinThetaO * cos2kAlpha[1] - cosThetaO * sin2kAlpha[1];
                cosThetaOp = cosThetaO * cos2kAlpha[1] + sinThetaO * sin2kAlpha[1];
            }

            // Handle remainder of p values for hair scale tilt
            else if (p == 1) {
                sinThetaOp = sinThetaO * cos2kAlpha[0] + cosThetaO * sin2kAlpha[0];
                cosThetaOp = cosThetaO * cos2kAlpha[0] - sinThetaO * sin2kAlpha[0];
            } else if (p == 2) {
                sinThetaOp = sinThetaO * cos2kAlpha[2] + cosThetaO * sin2kAlpha[2];
                cosThetaOp = cosThetaO * cos2kAlpha[2] - sinThetaO * sin2kAlpha[2];
            } else {
                sinThetaOp = sinThetaO;
                cosThetaOp = cosThetaO;
            }

            // Handle out-of-range cosThetaO from scale adjustment
            cosThetaOp = abs(cosThetaOp);
            pdf += Mp(cosThetaI, cosThetaOp, sinThetaI, sinThetaOp, v[p]) *
                apPdf[p] * Np(phi, p, s, gammaO, gammaT);
        }

        pdf += Mp(cosThetaI, cosThetaO, sinThetaI, sinThetaO, v[pMax]) *
            apPdf[pMax] * (1 / (2 * pi));

        return pdf;
    }


}  // namespace yocto

// -----------------------------------------------------------------------------
// IMPLEMENTATION OF SIGMA_A CONVERSION FUNCTIONS
// -----------------------------------------------------------------------------
namespace yocto {

    vec3f SigmaAFromConcentration(float ce, float cp){
        auto eumelaninSigmaA = vec3f{0.419f, 0.697f, 1.37f};
        auto pheomelaninSigmaA = vec3f{0.187f, 0.4f, 1.05f};

        auto sigma_a = ce * eumelaninSigmaA + cp * pheomelaninSigmaA;

        return sigma_a;
    }

    vec3f SigmaAFromReflectance(vec3f c, float beta_n){
        auto sigma_a = log(c) / 
                       (5.969f - 0.215f * beta_n + 2.532f * Sqr(beta_n) -
                        10.73f * Pow<3>(beta_n) + 5.574f *  Pow<4>(beta_n) +
                        0.245f *  Pow<5>(beta_n));

        sigma_a = sigma_a * sigma_a;

        return sigma_a;
    }

} // namespace yocto
