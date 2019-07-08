#pragma once

#include "Particle.hpp"


namespace picongpu
{
    namespace transitionRadiation
    {
        using complex_X = pmacc::math::Complex< float_X >;
        using complex_64 = pmacc::math::Complex< float_64 >;

        /** arbitrary margin which is necessary to prevent division by 0 error
         * created by particles moving in the plane of the foil.
         */
        float_X const DIV_BY_ZERO_MINIMUM = 1.e-7;

        class Calculator
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
            Calculator(
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
                    parMomPhi - detectorPhi, 
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

                float_X denominator = x * x - y * y;

                // Preventing division by 0
                if( math::abs( denominator ) < DIV_BY_ZERO_MINIMUM )
                {
                    if( denominator < 0.0 )
                        denominator = -DIV_BY_ZERO_MINIMUM;
                    else
                        denominator = DIV_BY_ZERO_MINIMUM;
                }

                return a * ( 1.0 / denominator );
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

                float_X denominator = x * x - y * y; 

                // Preventing division by 0
                if( math::abs( denominator ) < DIV_BY_ZERO_MINIMUM )
                {
                    if( denominator < 0.0 )
                        denominator = -DIV_BY_ZERO_MINIMUM;
                    else
                        denominator = DIV_BY_ZERO_MINIMUM;
                }

                return a * ( b - c ) * ( 1.0 / denominator );
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
                // If case for longitudinal moving particles... leads to 0 later in the kernel
                if ( math::abs( parMomCosTheta ) <= DIV_BY_ZERO_MINIMUM )
                    return complex_X( -1.0, 0.0 );
                    
                float_X const a = detectorSinTheta * parMomSinTheta * math::cos( parMomPhi - detectorPhi );
                float_X const b = - ( particle.getPosPara( ) ) * ( 1 / particle.getVel( ) - a / SPEED_OF_LIGHT) / ( parMomCosTheta );
                float_X const c = - detectorSinTheta * particle.getPosPerp( ) * math::cos( particle.getPosPhi( ) - detectorPhi );
                // float_X const a = 1.0;
                // float_X const b = 1.0;
                // float_X const c = 1.0;

                complex_X const fpara = complex_X( 0.0, b );
                complex_X const fperp = complex_X( 0.0, c );
                return fpara + fperp;
                
            }
        }; // class Calculator

        HDINLINE
        complex_X 
        formFactor(
            float_X const omega, 
            complex_X const exponent
        )
        {
            /* Does exactly what the name says */
            
            // If case for longitudinal moving particles
            if ( exponent.get_real() == -1.0 )
                return complex_X( 0.0, 0.0 );
            else
                return complex_X( 
                    math::exp( 
                        precisionCast< float_64 >( exponent * omega )
                    ) 
                );
        }
    } // namespace transitionRadiation
} // namespace picongpu