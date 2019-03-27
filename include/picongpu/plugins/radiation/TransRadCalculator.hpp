#pragma once

#include "CoolParticle.hpp"

namespace picongpu
{
    class TransRadCalculator
    {
        private:
        const CoolParticle& particle;
        const float3_X lookDirection;

        picongpu::float_X parMomSinTheta;
        picongpu::float_X parMomCosTheta;
        picongpu::float_X parMomSinPhi;
        picongpu::float_X parMomCosPhi;
        picongpu::float_X detectorSinTheta;
        picongpu::float_X detectorCosTheta;
        const picongpu::float_X parSqrtOnePlusUSquared;

        public: 
        HDINLINE
        TransRadCalculator(
            const CoolParticle& particle_set,
            const float3_X lookDirection
        ) : 
            particle(particle_set),
            lookDirection(lookDirection),
            parSqrtOnePlusUSquared( picongpu::math::sqrt( 1 + particle.getU() * particle.getU() ) )
        { 
            // Frequent calculations
            // Momentum Space for Particle:
            picongpu::math::sincos( particle.getMomTheta(), parMomSinTheta, parMomCosTheta );
            picongpu::math::sincos( particle.getMomPhi(), parMomSinPhi, parMomCosPhi );
            
            // Detector Position
            //const picongpu::float_X detectorTheta = picongpu::math::atan2(lookDirection.x(), lookDirection.y()) + picongpu::PI;
            const picongpu::float_X detectorTheta = picongpu::math::acos(lookDirection.y() / picongpu::math::sqrt((lookDirection * lookDirection).sumOfComponents()));
            picongpu::math::sincos( detectorTheta, detectorSinTheta, detectorCosTheta );
        }

        HDINLINE 
        picongpu::float_X calcEPerp(void) const
        {
            const picongpu::float_X uSquared = particle.getU() * particle.getU();
            const picongpu::float_X a = uSquared * parMomCosTheta * parMomSinTheta * parMomSinPhi * detectorCosTheta;
 
            // Denominator
            const picongpu::float_X x = parSqrtOnePlusUSquared - particle.getU() * parMomSinTheta * parMomCosPhi * detectorSinTheta;
            const picongpu::float_X y = particle.getU() * parMomCosTheta * detectorCosTheta;
            const picongpu::float_X xSquared = util::square<picongpu::float_X> ( x );
            const picongpu::float_X ySquared = util::square<picongpu::float_X> ( y );

            const picongpu::float_X denominator = xSquared - ySquared;

            return particle.getCharge() * a * (1.0 / denominator);
        }

        HDINLINE
        picongpu::float_X calcEPara(void) const
        {
            // const picongpu::float_X a = particle.getU() * parMomCosTheta;
            // const picongpu::float_X b = particle.getU() * parMomSinTheta * parMomCosPhi;
            // const picongpu::float_X c = parSqrtOnePlusUSquared * detectorSinTheta;

            // // Denominator
            // const picongpu::float_X x = parSqrtOnePlusUSquared - particle.getU() * parMomSinTheta * parMomCosPhi * detectorSinTheta;
            // const picongpu::float_X y = particle.getU() * parMomCosTheta * detectorCosTheta;
            // const picongpu::float_X xSquared = util::square<picongpu::float_X> ( x );
            // const picongpu::float_X ySquared = util::square<picongpu::float_X> ( y );

            // const picongpu::float_X denominator = xSquared - ySquared;
            // return particle.getCharge() * a * ( b - c ) * (1.0 / denominator);
            const picongpu::float_X a = particle.getU() * parMomCosTheta;
            const picongpu::float_X b = particle.getU() * parMomSinTheta * parMomCosPhi;
            const picongpu::float_X c = parSqrtOnePlusUSquared * detectorSinTheta;

            // Denominator
            const picongpu::float_X x = parSqrtOnePlusUSquared - particle.getU() * parMomSinTheta * parMomCosPhi * detectorSinTheta;
            const picongpu::float_X y = particle.getU() * parMomCosTheta * detectorCosTheta;
            const picongpu::float_X xSquared = util::square<picongpu::float_X> ( x );
            const picongpu::float_X ySquared = util::square<picongpu::float_X> ( y );

            const picongpu::float_X denominator = xSquared - ySquared;
            return particle.getCharge() * a * ( b - c ) * (1.0 / denominator);
        }
    };
}