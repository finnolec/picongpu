#pragma once

#include "CoolParticle.hpp"


namespace picongpu
{
    using complex_X = pmacc::math::Complex< float_X >;
    class TransRadCalculator
    {

        private:
        const CoolParticle& particle;
        const float3_X lookDirection;

        float_X parMomSinTheta;
        float_X parMomCosTheta;
        const float_X parMomPhi;
        float_X parMomSinPhi;
        float_X parMomCosPhi;
        float_X detectorSinTheta;
        float_X detectorCosTheta;
        const float_X detectorPhi;
        const float_X parSqrtOnePlusUSquared;

        public: 

        HDINLINE
        TransRadCalculator(
            const CoolParticle& particle_set,
            const float3_X lookDirection
        ) : 
            particle( particle_set ),
            lookDirection( lookDirection ),
            parMomPhi( particle.getMomPhi( ) ),
            detectorPhi( 
                picongpu::math::atan2(
                    lookDirection.z( ), 
                    lookDirection.x( )
                ) + picongpu::PI
            ),
            parSqrtOnePlusUSquared( 
                picongpu::math::sqrt( 1 + particle.getU( ) * particle.getU( ) ) 
            )
        { 
            // Frequent calculations
            // Momentum Space for Particle:
            picongpu::math::sincos( 
                particle.getMomTheta( ), 
                parMomSinTheta, 
                parMomCosTheta 
            );
            picongpu::math::sincos( 
                parMomPhi, 
                parMomSinPhi, 
                parMomCosPhi 
            );
            
            // Detector Position
            //const float_X detectorTheta = picongpu::math::atan2(lookDirection.x( ), lookDirection.y( )) + picongpu::PI;
            const float_X detectorTheta = picongpu::math::acos( 
                lookDirection.y( ) / picongpu::math::sqrt( ( 
                    lookDirection * lookDirection
                ).sumOfComponents( ) )
            );
            picongpu::math::sincos( 
                detectorTheta, 
                detectorSinTheta, 
                detectorCosTheta 
            );
        }

        HDINLINE 
        float_X calcEPerp( ) const
        {
            /* returns perpendicular part to movement direction of normalized energy
             * determined by formula:
             * E_perp = (u^2 cosPsi sinPsi sinPhi cosTheta) / 
             *          ((sqrt(1 + u^2) - u sinPsi cosPhi sinTheta)^2 - u^2 cosPsi^2 cosTheta^2)
             * where Psi is the azimuth angle of the particle momentum and theta is
             * the azimuth angle of the detector position to the movement direction y
             */
            const float_X uSquared = particle.getU( ) * particle.getU( );
            const float_X a = uSquared * parMomCosTheta * parMomSinTheta * parMomSinPhi * detectorCosTheta;
 
            // Denominator
            const float_X x = parSqrtOnePlusUSquared - particle.getU( ) * parMomSinTheta * parMomCosPhi * detectorSinTheta;
            const float_X y = particle.getU( ) * parMomCosTheta * detectorCosTheta;
            const float_X xSquared = util::square<float_X> ( x );
            const float_X ySquared = util::square<float_X> ( y );

            const float_X denominator = xSquared - ySquared;

            return particle.getCharge( ) * a * ( 1.0 / denominator );
        }

        HDINLINE
        float_X calcEPara( ) const
        {
            /* returns parallel part to movement direction of normalized energy
             * determined by formula:
             * E_perp = (u cosPsi (u sinPsi cosPhi - sqrt(1 + u^2) sinTheta)) / 
             *          ((sqrt(1 + u^2) - u sinPsi cosPhi sinTheta)^2 - u^2 cosPsi^2 cosTheta^2)
             * where Psi is the azimuth angle of the particle momentum and theta is
             * the azimuth angle of the detector position to the movement direction y
             */
            const float_X a = particle.getU( ) * parMomCosTheta;
            const float_X b = particle.getU( ) * parMomSinTheta * parMomCosPhi;
            const float_X c = parSqrtOnePlusUSquared * detectorSinTheta;

            // Denominator
            const float_X x = parSqrtOnePlusUSquared - particle.getU( ) * parMomSinTheta * parMomCosPhi * detectorSinTheta;
            const float_X y = particle.getU( ) * parMomCosTheta * detectorCosTheta;
            const float_X xSquared = util::square<float_X> ( x );
            const float_X ySquared = util::square<float_X> ( y );

            const float_X denominator = xSquared - ySquared;
            return particle.getCharge( ) * a * ( b - c ) * ( 1.0 / denominator );
        }

        HDINLINE
        complex_X calcFExponent( ) const
        {
            /* returns the exponent of the formfactor divided by \omega
             * this doesn't have a physical meaning, it's calculated here for performance reasons
             * determined by formula:
             * F_exp = - i z ( 1 / v - sinTheta sinPsi cos(Phi_P - Phi_D) / c ) / (cosPsi)
             *          - i sinTheta rho cos(Phi_P - Phi_D)
             */
            const float_X a = detectorSinTheta * parMomSinTheta * math::cos( parMomPhi - detectorPhi );
            const float_X b =  - (particle.getPosPara( ) - parameters::surfacePosition) * ( 1 / particle.getVel( ) - a / SPEED_OF_LIGHT) / ( parMomCosTheta);
            const float_X c = - detectorSinTheta * particle.getPosPerp( ) * math::cos( particle.getPosPhi( ) - detectorPhi);
            const complex_X fpara = complex_X( 0.0, b );
            const complex_X fperp = complex_X( 0.0, c );
            return fpara + fperp;
             
        }
    };

    HDINLINE
    complex_X makeExponentUseful(
        const float_X omega, 
        const  complex_X exponent
    )
    {
        /* Does exactly what the name says */
        return math::exp( exponent * omega );
    }
}