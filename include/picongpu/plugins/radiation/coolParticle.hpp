#pragma once

#include <pmacc/math/Vector.hpp>
#include "utilities.hpp"

namespace picongpu
{
    class CoolParticle
    {
        private:
        const vector_X momentum;
        const vector_X location;
        const picongpu::float_X charge;
        const picongpu::float_X mass;

        public:
        HDINLINE 
        CoolParticle( 
            const vector_X& location_set, 
            const vector_X& momentum_set,
            const picongpu::float_X charge_set,
            const picongpu::float_X mass_set
        ) :
            location(location_set),
            momentum(momentum_set),
            charge(charge_set),
            mass(mass_set)
        { }

        HDINLINE 
        vector_64 getMomentum(void) const
        {
            // returns momentum
            return momentum;
        }

        HDINLINE 
        picongpu::float_64 getU(void) const
        {
            // returns normalized momentum
            return calc_u(getMomentum());
        }

        HDINLINE 
        picongpu::float_64 getV(void) const
        {
            // return velocity
            return calc_velocity(getMomentum());
        }

        HDINLINE
        picongpu::float_64 getMomPhi(void) const
        {
            // returns polar angle phi of velocity
        }

        HDINLINE
        picongpu::float_64 getMomTheta(void) const
        {
            // returns azimuth angle theta of velocity
            return calcMomTheta(getMomentum());
        }

        private:
        HDINLINE 
        picongpu::float_64 calcU(const vector_X& momentum) const
        {
            //returns normalized momentum u = gama * beta
            const picongpu::float_32 gamma1 = calcGamma(momentum);
            const picongpu::float_32 beta1 = calcBeta(momentum);
            return gamma1 * beta1;
        }

        HDINLINE 
        picongpu::float_64 calcGamma(const vector_X& momentum) const
        {
            //returns lorentz factor gamma = sqrt(1/(1-beta**2))
            const picongpu::float_32 betaSquared = util::square<
                vector_X, 
                picongpu::float_32 
            > (calcBeta(momentum));
            return sqrt( 1.0 / ( 1 - betaSquared ));
        }

        HDINLINE 
        picongpu::float_64 calcBeta(const vector_X& momentum) const
        {
            // return beta = v / c
            return calcVelocity(momentum) * (1.0 / picongpu::SPEED_OF_LIGHT);
        }

        HDINLINE 
        picongpu::float_64 calcVelocity(const vector_X& momentum) const
        {
            //returns velocity v = p/m
            return std::sqrt(momentum * momentum) * (1.0 / mass);
        }

        HDINLINE
        picongpu::float_64 calcMomTheta(const vector_X& momentum) const
        {
            // returns azimuth angle between momentum and y-axis
            uY = calcU(momentum.y());
            uAbs = std::sqrt(calcU(momentum) * calcU(momentum));
            return std::acos(uY * (1.0 / uAbs));
        }

        HDINLINE
        picongpu::float_64 calcMomPhi(const vector_X& momentum) const
        {
            // returns polar angle on perpendicular plane
            uX = calcU(momentum.x());
            uY = calcU(momentum.y());
            return std::atan2(uY, uX) + picongpu::PI;
        }
    }; // class CoolParticle
} // namespace picongpu