#pragma once

#include "TransRadParticle.hpp"


namespace picongpu
{
    using complex_X = pmacc::math::Complex< float_X >;

    class TransRadCalculator
    {

    private:
        transitionRadiation::Particle const & particle;
        float3_X const lookDirection;

        float_X parMomSinTheta;
        float_X parMomCosTheta;
        float_X const parMomPhi;
        float_X parMomSinPhi;
        float_X parMomCosPhi;
        float_X detectorSinTheta;
        float_X detectorCosTheta;
        float_X const detectorPhi;
        float_X const parSqrtOnePlusUSquared;

    public: 

        HDINLINE
        TransRadCalculator(
            transitionRadiation::Particle const & particleSet,
            float3_X const lookDirection
        ) : 
            particle( particleSet ),
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
            float_X const detectorTheta = picongpu::math::acos( 
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
        float_X 
        calcEPerp( ) const
        {
            /* returns perpendicular part to movement direction of normalized energy
             * determined by formula:
             * E_perp = (u^2 cosPsi sinPsi sinPhi cosTheta) / 
             *          ((sqrt(1 + u^2) - u sinPsi cosPhi sinTheta)^2 - u^2 cosPsi^2 cosTheta^2)
             * where Psi is the azimuth angle of the particle momentum and theta is
             * the azimuth angle of the detector position to the movement direction y
             */
            float_X const uSquared = particle.getU( ) * particle.getU( );
            float_X const a = uSquared * parMomCosTheta * parMomSinTheta * 
                parMomSinPhi * detectorCosTheta;
 
            // Denominator
            float_X const x = parSqrtOnePlusUSquared - 
                particle.getU( ) * parMomSinTheta * parMomCosPhi * detectorSinTheta;
            float_X const y = particle.getU( ) * parMomCosTheta * detectorCosTheta;
            float_X const xSquared = util::square< float_X > ( x );
            float_X const ySquared = util::square< float_X > ( y );

            float_X const denominator = xSquared - ySquared;

            return particle.getCharge( ) * a * ( 1.0 / denominator );
        }

        HDINLINE
        float_X 
        calcEPara( ) const
        {
            /* returns parallel part to movement direction of normalized energy
             * determined by formula:
             * E_perp = (u cosPsi (u sinPsi cosPhi - sqrt(1 + u^2) sinTheta)) / 
             *          ((sqrt(1 + u^2) - u sinPsi cosPhi sinTheta)^2 - u^2 cosPsi^2 cosTheta^2)
             * where Psi is the azimuth angle of the particle momentum and theta is
             * the azimuth angle of the detector position to the movement direction y
             */
            float_X const a = particle.getU( ) * parMomCosTheta;
            float_X const b = particle.getU( ) * parMomSinTheta * parMomCosPhi;
            float_X const c = parSqrtOnePlusUSquared * detectorSinTheta;

            // Denominator
            float_X const x = parSqrtOnePlusUSquared - 
                particle.getU( ) * parMomSinTheta * parMomCosPhi * detectorSinTheta;
            float_X const y = particle.getU( ) * parMomCosTheta * detectorCosTheta;
            float_X const xSquared = util::square< float_X > ( x );
            float_X const ySquared = util::square< float_X > ( y );

            float_X const denominator = xSquared - ySquared;
            return particle.getCharge( ) * a * ( b - c ) * ( 1.0 / denominator );
        }

        HDINLINE
        complex_X 
        calcFExponent( ) const
        {
            /* returns the exponent of the formfactor divided by \omega
             * this doesn't have a physical meaning, it's calculated here for performance reasons
             * determined by formula:
             * F_exp = - i z ( 1 / v - sinTheta sinPsi cos(Phi_P - Phi_D) / c ) / (cosPsi)
             *          - i sinTheta rho cos(Phi_P - Phi_D)
             */
            float_X const a = detectorSinTheta * parMomSinTheta * math::cos( parMomPhi - detectorPhi );
            float_X const b =  - ( particle.getPosPara( ) - parameters::surfacePosition ) * ( 1 / particle.getVel( ) - a / SPEED_OF_LIGHT) / ( parMomCosTheta );
            float_X const c = - detectorSinTheta * particle.getPosPerp( ) * math::cos( particle.getPosPhi( ) - detectorPhi );
            complex_X const fpara = complex_X( 0.0, b );
            complex_X const fperp = complex_X( 0.0, c );
            return fpara + fperp;
             
        }
    };

    HDINLINE
    complex_X 
    makeExponentUseful(
        float_X const omega, 
        complex_X const exponent
    )
    {
        /* Does exactly what the name says */
        return math::exp( exponent * omega );
    }
}