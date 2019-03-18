#pragma once

#include <pmacc/math/Vector.hpp>
#include "utilities.hpp"

namespace picongpu
{
    class CoolParticle
    {
        public:
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
        { 
        }

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
            return calcU(getMomentum());
        }

        HDINLINE 
        picongpu::float_64 getV(void) const
        {
            // return velocity
            return calcVelocity(getMomentum());
        }

        HDINLINE 
        picongpu::float_64 calcGamma(const vector_X& momentum) const
        {
            //returns lorentz factor gamma = sqrt(1/(1-beta**2))
            const picongpu::float_64 betaSquared = util::square<
                picongpu::float_64, 
                picongpu::float_64 
            > (calcBeta(momentum));
            return picongpu::math::sqrt( 1.0 / ( 1 - betaSquared ));
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
        picongpu::float_64 calcU(const vector_X& momentum) const
        {
            //returns normalized momentum u = gama * beta
            const picongpu::float_64 gamma1 = calcGamma(momentum);
            const picongpu::float_64 beta1 = calcBeta(momentum);
            return gamma1 * beta1;
        }

        // Getters for Momentum in spherical coordinates
        HDINLINE
        picongpu::float_64 getMomPhi(void) const
        {
            //return polar angle phi of momentum
            return calcMomPhi();
        }

        HDINLINE
        picongpu::float_64 getMomTheta(void) const
        {
            //return azimuth angle psi of momentum
            return calcMomTheta();
        }

        HDINLINE
        picongpu::float_64 getMomAbs(void) const
        {
            //return absolute value of momentum
            return calcMomAbs();
        }

        private:
        // Calculators for Momentum in spherical coordinates
        HDINLINE
        picongpu::float_64 calcMomPhi(void) const
        {
            //return polar angle phi of momentum
            return picongpu::math::atan2(momentum.z(), momentum.x()) + picongpu::PI;
        }

        HDINLINE
        picongpu::float_64 calcMomTheta(void) const
        {
            //return azimuth angle psi of momentum
            return picongpu::math::acos(momentum.y() * (1.0 / getMomAbs()));
        }

        HDINLINE
        picongpu::float_64 calcMomAbs(void) const
        {
            //return absolute value of momentum
            return picongpu::math::sqrt(momentum * momentum);
        }
    }; // class CoolParticle
} // namespace picongpu