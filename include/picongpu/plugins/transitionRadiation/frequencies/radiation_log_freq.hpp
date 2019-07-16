/* Copyright 2013-2019 Heiko Burau, Rene Widera, Richard Pausch
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

#include "picongpu/simulation_defines.hpp"


namespace picongpu
{
    namespace transitionRadiation
    {
        namespace rad_log_frequencies
        {
            class FreqFunctor
            {
            public:
                FreqFunctor( void )
                {
                    omega_log_min = math::log( omega_min );
                    delta_omega_log = ( math::log( omega_max ) - omega_log_min ) / float_X( N_omega - 1 );
                }

                HDINLINE float_X operator( )( const int ID )
                {
                    return  math::exp( omega_log_min + ( float_X( ID ) ) * delta_omega_log );
                }

            private:
                float_X omega_log_min;
                float_X delta_omega_log;
            }; // FreqFunctor


            class InitFreqFunctor
            {
            public:
                InitFreqFunctor( void )
                { }

                HINLINE void Init( const std::string path )
                { }


                HINLINE FreqFunctor getFunctor( void )
                {
                    return FreqFunctor( );
                }
            }; // InitFreqFunctor


            /** Getter for parameters in output file
             * 
             * return parameters as string
             */
            HDINLINE
            std::string 
            getParameters( void )
            {
                std::string params = std::string( "log\t" );
                params += std::to_string( transitionRadiation::frequencies::N_omega ) + "\t";
                params += std::to_string( transitionRadiation::frequencies::SI::omega_min ) + "\t";
                params += std::to_string( transitionRadiation::frequencies::SI::omega_max ) + "\t";
                return params; 
            }
        } // namespace rad_log_frequencies
    } // namespace transitionRadiation
} // namespace picongpu
