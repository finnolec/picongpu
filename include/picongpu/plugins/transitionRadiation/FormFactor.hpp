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
#include "picongpu/plugins/radiation/utilities.hpp"


namespace picongpu
{
namespace plugins
{
namespace transitionRadiation
{
namespace macroParticleFormFactorbaseShape_3D
{
    /** general form factor class of discrete charge distribution of PIC particle shape of order T_shapeOrder
     *
     * @tparam T_shapeOrder order of charge distribution shape in PIC code used for radiation form factor
     */

    template< uint32_t T_shapeOrder >
    struct FormFactor
    {
        /** Form Factor for T_shapeOrder-order particle shape charge distribution of N discrete electrons:
         * \f[ | \mathcal{F} |^2 = N + (N*N - N) * (sinc^2(n_x * L_x * \omega) * sinc^2(n_y * L_y * \omega) * sinc^2(n_z * L_z * \omega))^T_shapeOrder \f]
         *
         * with observation direction (unit vector) \f$ \vec{n} = (n_x, n_y, n_z) \f$
         * and with:
         * @param N     = weighting
         * @param omega = frequency
         * @param L_d   = the size of the CIC-particle / cell in dimension d
         *
         * @param N = macro particle weighting
         * @param omega = frequency at which to calculate the  form factor
         * @param observer_unit_vec = observation direction
         * @return the Form Factor: \f$ \sqrt( | \mathcal{F} |^2 ) \f$
         */
        HDINLINE float_X operator()( const float_X N, const float_X omega, plugins::radiation::vector_X const & observer_unit_vec ) const
        {
            float_X sincValue = float_X( 1.0 );
            for( uint32_t d = 0; d < DIM3; ++d )
                sincValue *= math::sinc( observer_unit_vec[d] * cellSize[d] / ( SPEED_OF_LIGHT * float_X( 2.0 ) ) * omega );

            // here we combine sinc^2(..) with (...)^T_shapeOrder to ...^(2 * T_shapeOrder)
            return math::sqrt( N + ( N * N - N ) * plugins::radiation::util::details::pow( sincValue , 2 * T_shapeOrder ) );
        }
    };

} // namespace macroParticleFormFactorbaseShape_3D


namespace macroParticleFormFactorCIC_3D
{
    struct FormFactor : public macroParticleFormFactorbaseShape_3D::FormFactor< 1 >
    { };

} // namespace macroParticleFormFactorCIC_3D

namespace macroParticleFormFactorTSC_3D
{
    struct FormFactor : public macroParticleFormFactorbaseShape_3D::FormFactor< 2 >
    { };

} // namespace macroParticleFormFactorTSC_3D

namespace macroParticleFormFactorPCS_3D
{
    struct FormFactor : public macroParticleFormFactorbaseShape_3D::FormFactor< 3 >
    { };

} // namespace macroParticleFormFactorPCS_3D


namespace macroParticleFormFactorCIC_1Dy
{
    struct FormFactor
    {
        /** Form Factor for 1-d CIC charge distribution iy y of N discrete electrons:
         * \f[ | \mathcal{F} |^2 = N + (N*N - N) * sinc^2(n_y * L_y * \omega) \f]
         *
         * with observation direction (unit vector) \f$ \vec{n} = (n_x, n_y, n_z) \f$
         * and with:
         * @param N     = weighting
         * @param omega = frequency
         * @param L_d   = the size of the CIC-particle / cell in dimension d
         *
         * @param N = macro particle weighting
         * @param omega = frequency at which to calculate the  form factor
         * @param observer_unit_vec = observation direction
         * @return the Form Factor: \f$ \sqrt( | \mathcal{F} |^2 ) \f$
         */
        HDINLINE float_X operator()(const float_X N, const float_X omega, const plugins::radiation::vector_X observer_unit_vec) const
        {
            float_X const sinc = math::sinc( CELL_HEIGHT / ( SPEED_OF_LIGHT * float_X( 2.0 ) ) * omega );
            return math::sqrt( N + ( N * N - N ) * sinc * sinc );
        }
    };
    
} // macroParticleFormFactorCIC_1Dy


namespace macroParticleFormFactorGaussSpherical
{
    struct FormFactor
    {
        /** Form Factor for point-symmetric Gauss-shaped charge distribution of N discrete electrons:
        * \f[ <rho(r)> = N*q_e* 1/sqrt(2*pi*sigma^2) * exp(-0.5 * r^2/sigma^2) \f]
        * with sigma = 0.5*c/delta_t (0.5 because sigma is defined around center)
        *
        * @param N = macro particle weighting
        * @param omega = frequency at which to calculate the  form factor
        * @param observer_unit_vec = observation direction
        * @return the Form Factor: \f$ \sqrt( | \mathcal{F} |^2 ) \f$
        */
        HDINLINE float_X operator()(const float_X N, const float_X omega, const plugins::radiation::vector_X observer_unit_vec) const
        {
            /* currently a fixed sigma of DELTA_T * c is used to describe the distribution - might become a parameter */
            float_X const exponent = omega * float_X( 0.5 ) * DELTA_T;
            float_X const eFunction = math::exp( float_X( -0.5 ) * exponent * exponent );
            return math::sqrt( N + ( N * N - N ) * eFunction * eFunction );
        }
    };

} // macroParticleFormFactorGauss_spherical


namespace macroParticleFormFactorGaussCell
{
    struct FormFactor
    {
        /** Form Factor for per-dimension Gauss-shaped charge distribution of N discrete electrons:
        * \f[ <rho(r)> = N*q_e* product[d={x,y,z}](1/sqrt(2*pi*sigma_d^2) * exp(-0.5 * d^2/sigma_d^2)) \f]
        * with sigma_d = 0.5*cell_width_d*n_d
        *
        * @param N = macro particle weighting
        * @param omega = frequency at which to calculate the  form factor
        * @param observer_unit_vec = observation direction
        * @return the Form Factor: \f$ \sqrt( | \mathcal{F} |^2 ) \f$
        */
        HDINLINE float_X operator()(const float_X N, const float_X omega, const plugins::radiation::vector_X observer_unit_vec) const
        {
            float_X const widthExponent = observer_unit_vec.x() * CELL_WIDTH / ( SPEED_OF_LIGHT * float_X(2.0) ) * omega;
            float_X const heightExponent = observer_unit_vec.y() * CELL_HEIGHT / ( SPEED_OF_LIGHT * float_X(2.0) ) * omega;
            float_X const depthExponent = observer_unit_vec.z() * CELL_DEPTH / ( SPEED_OF_LIGHT * float_X(2.0) ) * omega;
            float_X const eFunction = math::exp( 
                float_X ( -0.5 ) * ( 
                    widthExponent * widthExponent +
                    heightExponent * heightExponent +
                    depthExponent * depthExponent
                )
            );
            return math::sqrt(
                N + ( N * N - N ) * eFunction * eFunction
            );
        }
    };

} // macroParticleFormFactorGauss_cell


namespace macroParticleFormFactorIncoherent
{
    struct FormFactor
    {
        /** Form Factor for an incoherent charge distribution:
        *
        * @param N = macro particle weighting
        * @param omega = frequency at which to calculate the  form factor
        * @param observer_unit_vec = observation direction
        * @return the Form Factor: \f$ \sqrt( | \mathcal{F} |^2 == \sqrt(weighting) \f$
        */
        HDINLINE float_X operator()(const float_X N, const float_X omega, const plugins::radiation::vector_X observer_unit_vec) const
        {
            return math::sqrt(N);
        }
    };

} // macroParticleFormFactorincoherent


namespace macroParticleFormFactorCoherent
{
    struct FormFactor
    {
        /** Form Factor for a coherent charge distribution:
        *
        * @param N = macro particle weighting
        * @param omega = frequency at which to calculate the  form factor
        * @param observer_unit_vec = observation direction
        * @return the Form Factor: \f$ \sqrt( | \mathcal{F} |^2 == \sqrt(weighting) \f$
        */
        HDINLINE float_X operator()(const float_X N, const float_X omega, const plugins::radiation::vector_X observer_unit_vec) const
        {
            return N;
        }
    };

} // macroParticleFormFactorcoherent
} // namespace transitionRadiation
} // namespace plugins
} // namespace picongpu
