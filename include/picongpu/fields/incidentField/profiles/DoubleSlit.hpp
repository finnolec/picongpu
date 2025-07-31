#pragma once

#include "picongpu/defines.hpp"
#include "picongpu/fields/incidentField/Functors.hpp"
#include "picongpu/fields/incidentField/profiles/PlaneWave.hpp"

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
                struct DoubleSlit
                {
                    //! Get text name of the incident field profile
                    HINLINE static std::string getName()
                    {
                        return "DoubleSlit";
                    }

                    static nlohmann::json metadata()
                    {
                        return T_Params::metadata();
                    }
                };

                namespace detail
                {
                    /** Unitless plane wave parameters
                     *
                     * @tparam T_Params user (SI) parameters
                     */
                    template<typename T_Params>
                    struct DoubleSlitUnitless : public PlaneWaveUnitless<T_Params>
                    {
                        //! User SI parameters
                        using Params = T_Params;

                        //! Base unitless parameters
                        using Base = PlaneWaveUnitless<T_Params>;

                        static constexpr float_X SLIT_DISTANCE
                            = static_cast<float_X>(Params::SLIT_DISTANCE_SI / sim.unit.length());

                        static constexpr float_X SLIT_WIDTH
                            = static_cast<float_X>(Params::SLIT_WIDTH_SI / sim.unit.length());
                    };

                    /** Plane wave incident E functor
                     *
                     * @tparam T_Params parameters
                     */
                    template<typename T_Params>
                    struct DoubleSlitFunctorIncidentE
                        : public DoubleSlitUnitless<T_Params>
                        , public PlaneWaveFunctorIncidentE<T_Params>
                    {
                    public:
                        //! Unitless parameters type
                        using Unitless = DoubleSlitUnitless<T_Params>;

                        //! Base class
                        using Base = PlaneWaveFunctorIncidentE<T_Params>;

                        /** Create a functor on the host side for the given time step
                         *
                         * @param currentStep current time step index, note that it is fractional
                         * @param unitField conversion factor from SI to internal units,
                         *                  fieldE_internal = fieldE_SI / unitField
                         */
                        HINLINE DoubleSlitFunctorIncidentE(float_X const currentStep, float3_64 const unitField)
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
                            return incidentField::detail::BaseSeparableFunctorE<T_Params>::operator()(
                                *this,
                                totalCellIdx);
                        }

                        /** Get position-dependent transversal scalar factor for the given position
                         *
                         * Interface required by Base.
                         *
                         * @param totalCellIdx cell index in the total domain (including all moving window slides)
                         */
                        HDINLINE float_X getTransversal(floatD_X const& totalCellIdx) const
                        {
                            float3_X pos = this->getInternalCoordinates(totalCellIdx);

                            // pos[0] is propagation direction
                            // pos[1] is polarization direction
                            auto const distanceLeft = pos[2] - Unitless::SLIT_DISTANCE / 2.0_X;
                            auto const distanceRight = pos[2] + Unitless::SLIT_DISTANCE / 2.0_X;

                            if(math::abs(distanceLeft) < Unitless::SLIT_WIDTH / 2.0_X
                               || math::abs(distanceRight) < Unitless::SLIT_WIDTH / 2.0_X)
                            // if (math::abs(pos[2]) < Unitless::SLIT_WIDTH / 2.0_X)
                            {
                                return 1.0_X;
                            }
                            else
                            {
                                return 0.0_X;
                            }
                        }
                    };
                } // namespace detail
            } // namespace profiles

            namespace traits::detail
            {
                /** Get type of incident field E functor for the plane wave profile type
                 *
                 * @tparam T_Params parameters
                 */
                template<typename T_Params>
                struct GetFunctorIncidentE<profiles::DoubleSlit<T_Params>>
                {
                    using type = profiles::detail::DoubleSlitFunctorIncidentE<T_Params>;
                };

                /** Get type of incident field B functor for the plane wave profile type
                 *
                 * @tparam T_Params parameters
                 */
                template<typename T_Params>
                struct GetFunctorIncidentB<profiles::DoubleSlit<T_Params>>
                {
                    using type = incidentField::detail::ApproximateIncidentB<
                        typename GetFunctorIncidentE<profiles::DoubleSlit<T_Params>>::type>;
                };
            } // namespace traits::detail
        } // namespace incidentField
    } // namespace fields
} // namespace picongpu
