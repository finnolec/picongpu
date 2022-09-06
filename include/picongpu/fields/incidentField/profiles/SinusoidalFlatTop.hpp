/* Copyright 2013-2022 Axel Huebl, Heiko Burau, Rene Widera, Richard Pausch, Sergei Bastrakov
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

#include "picongpu/fields/incidentField/Functors.hpp"
#include "picongpu/fields/incidentField/Traits.hpp"

#include <cstdint>
#include <string>


namespace picongpu
{
    namespace fields
    {
        namespace incidentField
        {
            namespace profiles
            {
                template<typename T_Params>
                struct SinusoidalFlatTop
                {
                    //! Get text name of the incident field profile
                    static HINLINE std::string getName()
                    {
                        return "SinusoidalFlatTop";
                    }
                };

                namespace detail
                {
                    /** Unitless plane wave parameters
                     *
                     * @tparam T_Params user (SI) parameters
                     */
                    template<typename T_Params>
                    struct SinusoidalFlatTopUnitless : public BaseParamUnitless<T_Params>
                    {
                        //! User SI parameters
                        using Params = T_Params;

                        //! Base unitless parameters
                        using Base = BaseParamUnitless<T_Params>;

                        // unit: UNIT_TIME
                        static constexpr float_X LASER_NOFOCUS_CONSTANT
                            = static_cast<float_X>(Params::LASER_NOFOCUS_CONSTANT_SI / UNIT_TIME);
                        // unit: UNIT_TIME
                        static constexpr float_X INIT_TIME = static_cast<float_X>(
                            (Params::RAMP_INIT * Params::PULSE_LENGTH_SI + Params::LASER_NOFOCUS_CONSTANT_SI)
                            / UNIT_TIME);
                        // unit: UNIT_LENGTH
                        static constexpr float_64 X_RAMP_LENGTH
                            = static_cast<float_X>(Params::X_RAMP_LENGTH_SI / UNIT_LENGTH);
                        static constexpr float_64 X_PLAT_MIN
                            = static_cast<float_X>(Params::X_PLAT_MIN_SI / UNIT_LENGTH);
                        static constexpr float_64 X_PLAT_MAX
                            = static_cast<float_X>(Params::X_PLAT_MAX_SI / UNIT_LENGTH);
                        // unit: UNIT_LENGTH
                        static constexpr float_64 Y_RAMP_LENGTH
                            = static_cast<float_X>(Params::Y_RAMP_LENGTH_SI / UNIT_LENGTH);
                        static constexpr float_64 Y_PLAT_MIN
                            = static_cast<float_X>(Params::Y_PLAT_MIN_SI / UNIT_LENGTH);
                        static constexpr float_64 Y_PLAT_MAX
                            = static_cast<float_X>(Params::Y_PLAT_MAX_SI / UNIT_LENGTH);
                        // unit: UNIT_LENGTH
                        static constexpr float_64 Z_RAMP_LENGTH
                            = static_cast<float_X>(Params::Z_RAMP_LENGTH_SI / UNIT_LENGTH);
                        static constexpr float_64 Z_PLAT_MIN
                            = static_cast<float_X>(Params::Z_PLAT_MIN_SI / UNIT_LENGTH);
                        static constexpr float_64 Z_PLAT_MAX
                            = static_cast<float_X>(Params::Z_PLAT_MAX_SI / UNIT_LENGTH);
                    };

                    /** Plane wave incident E functor
                     *
                     * @tparam T_Params parameters
                     */
                    template<typename T_Params>
                    struct SinusoidalFlatTopFunctorIncidentE
                        : public SinusoidalFlatTopUnitless<T_Params>
                        , public incidentField::detail::BaseSeparableFunctorE<T_Params>
                    {
                    public:
                        //! Unitless parameters type
                        using Unitless = SinusoidalFlatTopUnitless<T_Params>;

                        //! Base functor type
                        using Base = incidentField::detail::BaseSeparableFunctorE<T_Params>;

                        /** Create a functor on the host side for the given time step
                         *
                         * @param currentStep current time step index, note that it is fractional
                         * @param unitField conversion factor from SI to internal units,
                         *                  fieldE_internal = fieldE_SI / unitField
                         */
                        HINLINE SinusoidalFlatTopFunctorIncidentE(float_X const currentStep, float3_64 const unitField)
                            : Base(currentStep, unitField)
                        {
                        }

                        /** Calculate incident field E value for the given position
                         *
                         * Interface required by Base.
                         *
                         * @param totalCellIdx cell index in the total domain (including all moving window slides)
                         * @return incident field E value in internal units
                         */
                        HDINLINE float3_X operator()(floatD_X const& totalCellIdx) const
                        {
                            return Base::operator()(*this, totalCellIdx);
                        }

                        /** Get time-dependent longitudinal scalar factor for the given time
                         *
                         * @param time time moment to calculate the factor at
                         * @param phaseShift additional phase shift to add on top of everything else,
                         *                   in radian
                         */
                        HDINLINE float_X getLongitudinal(float_X const time, float_X const phaseShift) const
                        {
                            auto envelope = Unitless::AMPLITUDE;
                            auto const mue = 0.5_X * Unitless::RAMP_INIT * Unitless::PULSE_LENGTH;
                            auto const tau = Unitless::PULSE_LENGTH * math::sqrt(2.0_X);
                            auto const endUpramp = mue;
                            auto const startDownramp = mue + Unitless::LASER_NOFOCUS_CONSTANT;
                            auto integrationCorrectionFactor = 0.0_X;
                            if(time > startDownramp)
                            {
                                // downramp = end
                                auto const exponent = (time - startDownramp) / tau;
                                envelope *= exp(-0.5_X * exponent * exponent);
                                integrationCorrectionFactor = (time - startDownramp) / (Unitless::w * tau * tau);
                            }
                            else if(time < endUpramp)
                            {
                                // upramp = start
                                auto const exponent = (time - endUpramp) / tau;
                                envelope *= exp(-0.5_X * exponent * exponent);
                                integrationCorrectionFactor = (time - endUpramp) / (Unitless::w * tau * tau);
                            }

                            auto const timeOszi = time - endUpramp;
                            auto const phase = Unitless::w * timeOszi + Unitless::LASER_PHASE + phaseShift;
                            // to understand both components [sin(...) + t/tau^2 * cos(...)] see description above
                            return (math::sin(phase) + math::cos(phase) * integrationCorrectionFactor) * envelope;
                        }

                        /** Get position-dependent transversal scalar factor for the given position
                         *
                         * Interface required by Base.
                         *
                         * @param totalCellIdx cell index in the total domain (including all moving window slides)
                         */
                        HDINLINE float_X getTransversal(floatD_X const& totalCellIdx) const
                        { 
                            float3_X const pos = this->getInternalCoordinates(totalCellIdx);
                            float_X transversalFactor = 1.0;

                            if((Unitless::DIR_X != 1.0)){
                                if((pos[0] <= (Unitless::X_PLAT_MIN - Unitless::X_RAMP_LENGTH)) || (pos[0] >= (Unitless::X_PLAT_MAX + Unitless::X_RAMP_LENGTH))){
                                    transversalFactor = 0.0;
                                } else if ((pos[0] > (Unitless::X_PLAT_MIN - Unitless::X_RAMP_LENGTH)) && (pos[0] < Unitless::X_PLAT_MIN)) {
                                    transversalFactor *= 0.5; //(1.0 - math::cos(
                                        //PI * (pos[0] - Unitless::X_PLAT_MIN + Unitless::X_RAMP_LENGTH) / Unitless::X_RAMP_LENGTH)) / 2.0;
                                } else if ((pos[0] > Unitless::X_PLAT_MAX) && (pos[0] < (Unitless::X_PLAT_MAX + Unitless::X_RAMP_LENGTH))) {
                                    transversalFactor *= 0.5; //(1.0 + math::cos(PI * (pos[0] - Unitless::X_PLAT_MAX) / Unitless::X_RAMP_LENGTH)) / 2.0;
                                }
                            }
/*
                            if((Unitless::DIR_Y != 1.0) && (transversalFactor != 0.0)){
                                if((pos[1] <= (Unitless::Y_PLAT_MIN - Unitless::Y_RAMP_LENGTH)) || (pos[1] >= (Unitless::Y_PLAT_MAX + Unitless::Y_RAMP_LENGTH))){
                                    transversalFactor = 0.0;
                                } else if ((pos[1] > (Unitless::Y_PLAT_MIN - Unitless::Y_RAMP_LENGTH)) && (pos[1] < Unitless::Y_PLAT_MIN)) {
                                    transversalFactor *= 0.5; //(1.0 - math::cos(
                                        //PI * (pos[1] - Unitless::Y_PLAT_MIN + Unitless::Y_RAMP_LENGTH) / Unitless::Y_RAMP_LENGTH)) / 2.0;
                                } else if ((pos[1] > Unitless::Y_PLAT_MAX) && (pos[1] < (Unitless::Y_PLAT_MAX + Unitless::Y_RAMP_LENGTH))) {
                                    transversalFactor *= 0.5; //(1.0 + math::cos(PI * (pos[1] - Unitless::Y_PLAT_MAX) / Unitless::Y_RAMP_LENGTH)) / 2.0;
                                }
                            }

                            if((Unitless::DIR_Z != 1.0) && (transversalFactor != 0.0)){
                                if((pos[2] <= (Unitless::Z_PLAT_MIN - Unitless::Z_RAMP_LENGTH)) || (pos[2] >= (Unitless::Z_PLAT_MAX + Unitless::Z_RAMP_LENGTH))){
                                    transversalFactor = 0.0;
                                } else if ((pos[2] > (Unitless::Z_PLAT_MIN - Unitless::Z_RAMP_LENGTH)) && (pos[2] < Unitless::Z_PLAT_MIN)) {
                                    transversalFactor *= 0.5; //(1.0 - math::cos(
                                        //PI * (pos[2] - Unitless::Z_PLAT_MIN + Unitless::Z_RAMP_LENGTH) / Unitless::Z_RAMP_LENGTH)) / 2.0;
                                } else if ((pos[2] > Unitless::Z_PLAT_MAX) && (pos[2] < (Unitless::Z_PLAT_MAX + Unitless::Z_RAMP_LENGTH))) {
                                    transversalFactor *= 0.5; //(1.0 + math::cos(PI * (pos[2] - Unitless::Z_PLAT_MAX) / Unitless::Z_RAMP_LENGTH)) / 2.0;
                                }
                            }
*/
                            return transversalFactor;
                        }
                    };
                } // namespace detail
            } // namespace profiles

            namespace detail
            {
                /** Get type of incident field E functor for the plane wave profile type
                 *
                 * @tparam T_Params parameters
                 */
                template<typename T_Params>
                struct GetFunctorIncidentE<profiles::SinusoidalFlatTop<T_Params>>
                {
                    using type = profiles::detail::SinusoidalFlatTopFunctorIncidentE<T_Params>;
                };

                /** Get type of incident field B functor for the plane wave profile type
                 *
                 * @tparam T_Params parameters
                 */
                template<typename T_Params>
                struct GetFunctorIncidentB<profiles::SinusoidalFlatTop<T_Params>>
                {
                    using type = detail::ApproximateIncidentB<
                        typename GetFunctorIncidentE<profiles::SinusoidalFlatTop<T_Params>>::type>;
                };
            } // namespace detail
        } // namespace incidentField
    } // namespace fields
} // namespace picongpu
