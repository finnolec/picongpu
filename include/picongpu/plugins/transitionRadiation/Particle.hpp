#pragma once

namespace picongpu
{
    namespace transitionRadiation
    {
        /** Particle class for transition radiation calculation.
         * 
         * @param locationSet global position of the macro-particle
         * @param momentumSet momentum of macro-particle
         * @param charge
         */
        class Particle
        {
        private:
            float3_X const & momentum;
            float3_X const & location;
            float_X const mass;

        public:
            HDINLINE 
            Particle( 
                float3_X const & locationSet, 
                float3_X const & momentumSet,
                float_X const massSet
            ) :
                location( locationSet ),
                momentum( momentumSet ),
                mass( massSet )
            { 
            }

            //! @return momentum
            HDINLINE 
            float3_X 
            getMomentum( ) const
            {
                return momentum;
            }

            //! @return normalized momentum
            HDINLINE 
            float_X 
            getU( ) const
            {
                return calcU( getMomentum( ) );
            }

            //! @return velocity v = beta * c
            HDINLINE
            float_X 
            getVel( ) const
            {
                return calcBeta( getMomentum( ) ) * picongpu::SPEED_OF_LIGHT;
            }

            //! @return polar angle phi of momentum
            HDINLINE
            float_X 
            getMomPhi( ) const
            {
                return calcMomPhi( );
            }

            //! @return azimuth angle psi of momentum
            HDINLINE
            float_X 
            getMomTheta( ) const
            {
                return calcMomTheta( );
            }

            //! @return absolute value of momentum
            HDINLINE
            float_X 
            getMomAbs( ) const
            {
                return calcMomAbs( );
            }
            
            //! @return radial, perpendicular component of location in cylindrical coordinates
            HDINLINE
            float_X 
            getPosPerp( ) const
            {
                return calcPosRho( );
            }

            //! @return parallel component to z of location in cylindrical coordinates
            HDINLINE
            float_X 
            getPosPara( ) const
            {
                return location.y( );
            }

            //! @return polar angle of location in cylindrical coordinates
            HDINLINE
            float_X 
            getPosPhi( ) const
            {
                return calcPosPhi( );
            }

        private:
            //! @return beta=v/c
            HDINLINE 
            float_X 
            calcBeta(
                float3_X const & momentum
            ) const
            {
                float_X const gamma = calcGamma(momentum);
                return picongpu::math::sqrt(1 - 1 / (gamma * gamma) );
            }

            //! @return gamma = E/(mc^2)
            HDINLINE 
            float_X 
            calcGamma(
                float3_X const & momentum
            ) const
            {
                float_X const x = 1.0 / ( mass * picongpu::SPEED_OF_LIGHT ) * ( mass * picongpu::SPEED_OF_LIGHT );
                return picongpu::math::sqrt( 
                    1.0 + ( momentum * momentum ).sumOfComponents( ) * x 
                );
            }

            //! @return normalized momentum u = gamma * beta
            HDINLINE 
            float_X 
            calcU(
                float3_X const & momentum
            ) const
            {
                float_X const gamma1 = calcGamma( momentum );
                float_X const beta1 = calcBeta( momentum );
                return gamma1 * beta1;
            }

            //! @return polar angle phi of momentum
            HDINLINE
            float_X 
            calcMomPhi( ) const
            {
                return picongpu::math::atan2(
                    momentum.x( ), 
                    momentum.z( )
                ) + picongpu::PI;
            }

            //! @return azimuth angle psi of momentum
            HDINLINE
            float_X 
            calcMomTheta( ) const
            {
                //because of floating point precision x^2+y^2+z^2<y^2 for x,z<<z
                float_X const momAbs = getMomAbs( );
                if( momAbs <= momentum.y( ) )
                    return 0.0;
                else
                    return picongpu::math::acos( momentum.y( ) * ( 1.0 / momAbs ) );
            }

            //! @return absolute value of momentum
            HDINLINE
            float_X 
            calcMomAbs( ) const
            {
                return picongpu::math::sqrt( 
                    momentum.x( ) * momentum.x( ) + 
                    momentum.y( ) * momentum.y( ) + 
                    momentum.z( ) * momentum.z( ) 
                );
            }

            //! @return radial component of location in cylindrical coordinates
            HDINLINE
            float_X 
            calcPosRho( ) const
            {
                return picongpu::math::sqrt(
                    location.x( ) * location.x( ) + 
                    location.z( ) * location.z( )
                );
            }

            //! @return radial angle phi of location in cylindrical coordinates
            HDINLINE
            float_X 
            calcPosPhi( ) const
            {
                return picongpu::math::atan2(
                    location.x( ), 
                    location.z( )
                ) + picongpu::PI;
            }
        }; // class Particle
    } // namespace transitionRadiation
} // namespace picongpu