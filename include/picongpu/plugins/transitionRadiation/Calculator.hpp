#pragma once

#include "Particle.hpp"


namespace picongpu
{
    namespace transitionRadiation
    {
        using complex_X = pmacc::math::Complex< float_X >;
        using complex_64 = pmacc::math::Complex< float_64 >;

        /* Arbitrary margin which is necessary to prevent division by 0 error
         * created by particles moving in the plane of the foil.
         */
        float_X const DIV_BY_ZERO_MINIMUM = 1.e-7;

        /** Calculator class for calculation of transition radiation.
         * 
         * @param particleSet transitionRadiation::Particle to compute transition radiation for
         * @param lookDirection vector with observation direction
         */
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

            /** Perpendicular part of normalized energy
             * 
             * Calculates perpendicular part to movement direction of normalized energy
             * determined by formula:
             * @f[E_{perp} = (u^2 \cos{\psi} \sin{\psi} \sin{\phi} \cos{\theta}) / 
             *          ((\sqrt{1 + u^2} - u \sin{\psi} \cos{\phi} \sin{\theta})^2 - u^2 \cos{\phi}^2 \cos{\theta}^2)@f]
             * where Psi is the azimuth angle of the particle momentum and theta is
             * the azimuth angle of the detector position to the movement direction y
             * 
             * @return perpendicular part of normalized energy
             */
            HDINLINE 
            float_X 
            calcEPerp( ) const
            {
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

            /** Parallel part of normalized energy
             * 
             * Calculates parallel part to movement direction of normalized energy
             * determined by formula:
             * @f[E_{para} = (u \cos{\psi} (u \sin{\psi} \cos{\phi} - \sqrt{1 + u^2} \sin{\theta}) / 
             *          ((\sqrt{1 + u^2} - u \sin{\psi} \cos{\phi} \sin{\theta})^2 - u^2 \cos{\phi}^2 \cos{\theta}^2)@f]
             * where Psi is the azimuth angle of the particle momentum and theta is
             * the azimuth angle of the detector position to the movement direction y
             * 
             * @return parallel part of normalized energy
             */
            HDINLINE
            float_X 
            calcEPara( ) const
            {
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

            /** Exponent of form factor
             * 
             * Calculates the exponent of the formfactor divided by \omega
             * It represents the phase of a single electron in the bunch, but it is mostly
             * calculated for performance reasons.
             * \f[ F_exp = - i z ( 1 / v - \sin{\theta} \sin{\psi} \cos{\phi_P - \phi_D} / c ) / \cos{\phi}
             *          - i \sin{\theta} \rho \cos{\phi_P - \phi_D} \f]
             * 
             */
            HDINLINE
            complex_X 
            calcFExponent( ) const
            {
                // If case for longitudinal moving particles... leads to 0 later in the kernel
                if ( math::abs( parMomCosTheta ) <= DIV_BY_ZERO_MINIMUM )
                    return complex_X( -1.0, 0.0 );
                    
                float_X const a = detectorSinTheta * parMomSinTheta * math::cos( parMomPhi - detectorPhi );
                float_X const b = - ( particle.getPosPara( ) ) * ( 1 / particle.getVel( ) - a / SPEED_OF_LIGHT) / ( parMomCosTheta );
                float_X const c = - detectorSinTheta * particle.getPosPerp( ) * math::cos( particle.getPosPhi( ) - detectorPhi );

                complex_X const fpara = complex_X( 0.0, b );
                complex_X const fperp = complex_X( 0.0, c );
                return fpara + fperp;
                
            }
        }; // class Calculator

        /** Formfactor
         * 
         * Calculates of the electron bunch with the exponent calculated by the 
         * Calculator class. 
         * 
         * @f[F = \exp{ F_{exp} * \omega }@f]
         * 
         * @param omega observed frequency
         * @param exponent exponent of exponential function
         */
        HDINLINE
        complex_X 
        formFactor(
            float_X const omega, 
            complex_X const exponent
        )
        {
            // If case for longitudinal moving particles
            if ( exponent.get_real() == -1.0 )
                return complex_X( 0.0, 0.0 );
            else
                return complex_X( 
                    math::exp( 
                        //precisionCast< float_64 >( exponent * omega )
                        exponent * omega
                    ) 
                );
        }
    } // namespace transitionRadiation
} // namespace picongpu