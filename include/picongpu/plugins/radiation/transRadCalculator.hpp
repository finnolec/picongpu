#pragma once

#include "coolParticle.hpp"

namespace picongpu
{
    class TransRadCalculator
    {
        private:
        const CoolParticle& particle;
        const vector_64 lookDirection;

        public: 
        HDINLINE
        TransRadCalculator(
            const CoolParticle& particle_set,
            const vector_64 lookDirection
        ) : 
            particle(particle_set),
            lookDirection(lookDirection)
        { 
        }

        HDINLINE 
        picongpu::float_64 calcEPerp(void) const
        {
            // sine and cosine calculations for spherical momentum coordinates from Particle
            picongpu::float_64 parMomCosTheta;
            picongpu::float_64 parMomSinTheta;
            picongpu::math::sincos( particle.getMomTheta(), parMomSinTheta, parMomCosTheta );
            
            picongpu::float_64 parMomSinPhi;
            picongpu::float_64 parMomCosPhi;
            picongpu::math::sincos(particle.getMomPhi(), parMomSinPhi, parMomCosPhi);

            const picongpu::float_64 parSqrtOnePlusUSquared = picongpu::math::sqrt( 1 + util::square<picongpu::float_32>(particle.getU()) );

            // sine and cosine for detector position
            const picongpu::float_64 detectorTheta = picongpu::math::atan2( lookDirection.y(), lookDirection.x() ) + picongpu::PI;
            picongpu::float_64 detectorSinTheta;
            picongpu::float_64 detectorCosTheta;
            picongpu::math::sincos(detectorTheta, detectorSinTheta, detectorCosTheta);

            const picongpu::float_64 uSquared = util::square<picongpu::float_64> ( particle.getU() );

            // Denominator
            const picongpu::float_64 x = parSqrtOnePlusUSquared - particle.getU() * parMomSinPhi * parMomCosPhi * detectorSinTheta;
            const picongpu::float_64 y = particle.getU() * parMomSinTheta * detectorCosTheta;
            const picongpu::float_64 xSquared = util::square<picongpu::float_64> ( x );

            const picongpu::float_64 denominator = xSquared - y;

            return uSquared * parMomCosTheta * parMomSinTheta * parMomSinPhi * detectorCosTheta * (1.0 / denominator);
        }

        HDINLINE
        picongpu::float_64 calcEPara(void) const
        {
            // sine and cosine calculations for spherical momentum coordinates from Particle
            picongpu::float_64 parMomSinTheta;
            picongpu::float_64 parMomCosTheta;
            picongpu::math::sincos( particle.getMomTheta(), parMomSinTheta, parMomCosTheta );

            picongpu::float_64 parMomSinPhi;
            picongpu::float_64 parMomCosPhi;
            picongpu::math::sincos( particle.getMomPhi(), parMomSinPhi, parMomCosPhi );

            const picongpu::float_64 parSqrtOnePlusUSquared = picongpu::math::sqrt( 1 + util::square<picongpu::float_32>(particle.getU()) );
        
            // sine and cosine for detector position
            const picongpu::float_64 detectorTheta = picongpu::math::atan2( lookDirection.y(), lookDirection.x() ) + picongpu::PI;
            picongpu::float_64 detectorSinTheta;
            picongpu::float_64 detectorCosTheta;
            picongpu::math::sincos( detectorTheta, detectorSinTheta, detectorCosTheta );

            const picongpu::float_64 a = particle.getU() * parMomCosTheta;
            const picongpu::float_64 b = particle.getU() * parMomSinTheta * parMomCosPhi;
            const picongpu::float_64 c = picongpu::math::sqrt( 1 + util::square<picongpu::float_32>(particle.getU()) ) * detectorSinTheta;

            // Denominator
            const picongpu::float_64 x = parSqrtOnePlusUSquared - particle.getU() * parMomSinPhi * parMomCosPhi * detectorSinTheta;
            const picongpu::float_64 y = particle.getU() * parMomSinTheta * detectorCosTheta;
            const picongpu::float_64 xSquared = util::square<picongpu::float_64> ( x );

            const picongpu::float_64 denominator = xSquared - y;

            return a * ( b - c ) * (1.0 / denominator);
        }
    };
}