#pragma once

#include "utilities.hpp"


namespace picongpu
{
    class CoolParticle
    {
        public:
        const float3_X& momentum;
        const float3_X& location;
        const float_X charge;
        const float_X mass;

        public:
        HDINLINE 
        CoolParticle( 
            const float3_X& location_set, 
            const float3_X& momentum_set,
            const float_X charge_set,
            const float_X mass_set
        ) :
            location( location_set ),
            momentum( momentum_set ),
            charge( charge_set ),
            mass( mass_set )
        { 
        }

        HDINLINE 
        float3_X getMomentum( ) const
        {
            // returns momentum
            return momentum;
        }

        HDINLINE
        float_X getCharge( ) const
        {
            return charge;
        }

        HDINLINE 
        float_X getU( ) const
        {
            // returns normalized momentum
            return calcU( getMomentum( ) );
        }

        HDINLINE
        float_X getVel( ) const
        {
            return calcBeta( getMomentum( ) ) * picongpu::SPEED_OF_LIGHT;
        }

        // Getters for Momentum in spherical coordinates
        HDINLINE
        float_X getMomPhi( ) const
        {
            //return polar angle phi of momentum
            return calcMomPhi( );
        }

        HDINLINE
        float_X getMomTheta( ) const
        {
            //return azimuth angle psi of momentum
            return calcMomTheta( );
        }

        HDINLINE
        float_X getMomAbs( ) const
        {
            //return absolute value of momentum
            return calcMomAbs( );
        }

        HDINLINE
        float_X getPosPerp( ) const
        {
            // return radial, perpendicular component of location in cylindrical coordinates
            return calcPosRho( );
        }

        HDINLINE
        float_X getPosPara( ) const
        {
            // return parallel component to z of location in cylindrical coordinates
            return location.y( );
        }

        HDINLINE
        float_X getPosPhi( ) const
        {
            //return polar angle of location in cylindrical coordinates
            return calcPosPhi( );
        }

        private:
        // Calculators for Momentum in spherical coordinates
        HDINLINE float_X calcBeta(
            const float3_X& momentum
        ) const
        {
            // returns beta=v/c
            const float_X gamma1 = calcGamma( momentum );
            const float_X x = util::square<float_X, float_X >(
                ( 1.0 / ( mass * picongpu::SPEED_OF_LIGHT * gamma1 ) )
            );
            return picongpu::math::sqrt( ( momentum * momentum ).sumOfComponents( ) * x);
        }

        HDINLINE float_X calcGamma(
            const float3_X& momentum
        ) const
        {
            // return gamma = E/(mc^2)
            const float_X x = util::square< float_X, float_X > ( ( 
                1.0 / ( mass * picongpu::SPEED_OF_LIGHT ) ) 
            );
            return picongpu::math::sqrt( 
                1.0 + ( momentum * momentum ).sumOfComponents( ) * x 
            );
        }

        HDINLINE 
        float_X calcU(
            const float3_X& momentum
        ) const
        {
            //returns normalized momentum u = gama * beta
            const float_X gamma1 = calcGamma( momentum );
            const float_X beta1 = calcBeta( momentum );
            return gamma1 * beta1;
        }

        HDINLINE
        float_X calcMomPhi( ) const
        {
            //return polar angle phi of momentum
            return picongpu::math::atan2(
                momentum.z( ), 
                momentum.x( )
            ) + picongpu::PI;
        }

        HDINLINE
        float_X calcMomTheta( ) const
        {
            //return azimuth angle psi of momentum
            //because of floating point precision x^2+y^2+z^2<y^2 for x,z<<z
            const float_X momAbs = getMomAbs( );
            if( momAbs <= momentum.y( ) )
            {
                return 0.0;
            }
            else
            {
                return picongpu::math::acos( momentum.y( ) * ( 1.0 /momAbs ) );
            }
        }

        HDINLINE
        float_X calcMomAbs( ) const
        {
            //return absolute value of momentum
            return picongpu::math::sqrt( 
                momentum.x( ) * momentum.x( ) + 
                momentum.y( ) * momentum.y( ) + 
                momentum.z( ) * momentum.z( ) 
            );
        }

        HDINLINE
        float_X calcPosRho( ) const
        {
            // return radial component of location in cylindrical coordinates
            return picongpu::math::sqrt(
                location.x( ) * location.x( ) + 
                location.z( ) * location.z( )
            );
        }

        HDINLINE
        float_X calcPosPhi( ) const
        {
            // return radial angle phi of location in cylindrical coordinates
            return picongpu::math::atan2(
                location.z( ), 
                location.x( )
            ) + picongpu::PI;
        }

    }; // class CoolParticle
} // namespace picongpu