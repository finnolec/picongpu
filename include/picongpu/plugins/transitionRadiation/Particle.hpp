/* Copyright 2013-2019 Heiko Burau, Rene Widera, Richard Pausch, Finn-Ole Carstens
 *
 * This file is part of PIConGPU.
 *
 * PIConGPU is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * PIConGPU is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with PIConGPU.
 * If not, see <http://www.gnu.org/licenses/>.
 */

#pragma once

namespace picongpu
{
namespace plugins
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
            float_X const gamma = calcGamma( momentum );
            float_X const beta = calcBeta( momentum );
            return gamma * beta;
        }

        //! @return velocity v = beta * c
        HDINLINE
        float_X
        getVel( ) const
        {
            return calcBeta( momentum ) * picongpu::SPEED_OF_LIGHT;
        }

        //! @return polar angle phi of momentum
        HDINLINE
        float_X
        getMomPhi( ) const
        {
            // add pi to atan2 function, because phi is in range from 0 to 2 pi
            return picongpu::math::atan2(
                momentum.x( ),
                momentum.z( )
            ) + picongpu::PI;
        }

        //! @return azimuth angle psi of momentum
        HDINLINE
        float_X
        getMomTheta( ) const
        {
            //because of floating point precision x^2+y^2+z^2<y^2 for x,z<<y
            float_X const momAbs = getMomAbs( );
            if( momAbs <= momentum.y( ) )
                return 0.0;
            else
                return picongpu::math::acos( momentum.y( ) * ( 1.0 / momAbs ) );
        }

        //! @return absolute value of momentum
        HDINLINE
        float_X
        getMomAbs( ) const
        {
            return picongpu::math::sqrt(
                momentum.x( ) * momentum.x( ) +
                momentum.y( ) * momentum.y( ) +
                momentum.z( ) * momentum.z( )
            );
        }

        //! @return radial, perpendicular component of location in cylindrical coordinates
        HDINLINE
        float_X
        getPosPerp( ) const
        {
            return picongpu::math::sqrt(
                location.x( ) * location.x( ) +
                location.z( ) * location.z( )
            );
        }

        //! @return parallel component to y of location in cylindrical coordinates
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
            // add pi to atan2 function, because phi is in range from 0 to 2 pi
            return picongpu::math::atan2(
                location.x( ),
                location.z( )
            ) + picongpu::PI;
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
            return picongpu::math::sqrt(1.0 - 1.0 / (gamma * gamma) );
        }

        //! @return gamma = E/(mc^2)
        HDINLINE
        float_X
        calcGamma(
            float3_X const & momentum
        ) const
        {
            float_X const massTimesC = mass * picongpu::SPEED_OF_LIGHT;
            return picongpu::math::sqrt(
                1.0 + ( momentum * momentum ).sumOfComponents( ) / ( massTimesC * massTimesC )
            );
        }
    }; // class Particle

} // namespace transitionRadiation
} // namespace plugins
} // namespace picongpu
