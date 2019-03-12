#pragma once

#include "coolParticle.hpp"

namespace picongpu
{
    class TransRadCalculator
    {
        private:
        const Particle& particle;
        const picongpu::float_64 detectorPhi;
        const picongpu::float_64 detectorTheta;

        const picongpu::float_64 parMomSinTheta;
        const picongpu::float_64 parMomCosTheta;
        const picongpu::float_64 parMomSinPhi;
        const picongpu::float_64 parMomCosPhi;
        const picongpu::float_64 parSqrtOnePlusUSquared;

        const picongpu::float_64 detSinTheta;
        const picongpu::float_64 detCosTheta;

        public: 
        HDINLINE
        TransRadCalculator(
            const CoolParticcle& particle_set,
            const vector_64 lookDirection
        ) : 
            particle(particle_set)
        { 
            const lookPhi = std::acos( lookDirection.y() )
            const lookTheta = atan2( lookDirection.y(), lookDirection.x() ) + picongpu::PI;
            const parMomSinTheta = std::sin( particle.getMomTheta() );
            const parMomCosTheta = std::cos( particle.getMomTheta() );
            const parMomSinPhi = std::sin( particle.getMomPhi() );
            const parMomCosPhi = std::cos( particle.getMomPhi() );
            const parSqrtOnePlusUSquared = std::sqrt( 1 + util::square<picongpu::float_32>(particle.getU()) );
            const lookSinTheta = std::sin( lookTheta );
            const lookCosTheta = std::cos( lookTheta );
        }

        HDINLINE 
        picongpu::float_64 calcEPerp(void)
        {
            const picongpu::float_64 uSquared = util::square<picongpu::float_64> ( particle.getU() );
            return uSquared * parMomCosTheta * parMomSinTheta * parMomSinPhi * lookCosTheta * (1.0 / calcDenominator);
        }

        HDINLINE
        picongpu::float_64 calcEPara(void)
        {
            const picongpu::float_64 a = particle.getU() * parMomCosTheta;
            const picongpu::float_64 b = particle.getU() * parMomSinTheta * parMomCosPhi;
            const picongpu::float_64 c = parSqrtOnePlusUSquared * detSinTheta;
            return a * ( b - c ) / calcDenominator;
        }

        private:

        HDINLINE
        picongpu::float_64 calcDenominator(void)
        {
            const picongpu::float_64 x = parSqrtOnePlusUSquared - particle.getU() * parMomSinPhi * parMomCosPhi * lookSinTheta;
            const picongpu::float_64 y = particle.getU() * parMomSinTheta * lookCosTheta;
            const picongpu::float_64 xSquared = util::square<picongpu::float_64> ( x );
            return xSquared - y;
        }
    }
}