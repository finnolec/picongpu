#pragma once

#include "utilities.hpp"

namespace picongpu
{
    class CoolParticle
    {
        public:
        const float3_X& momentum;
        const float3_X& location;
        const picongpu::float_X charge;
        const picongpu::float_X mass;

        public:
        HDINLINE 
        CoolParticle( 
            const float3_X& location_set, 
            const float3_X& momentum_set,
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
        float3_X getMomentum(void) const
        {
            // returns momentum
            return momentum;
        }

        HDINLINE 
        picongpu::float_X getU(void) const
        {
            // returns normalized momentum
            return calcU(getMomentum());
        }

        HDINLINE
        picongpu::float_X getCharge(void) const
        {
            return charge;
        }
/*
        HDINLINE 
        picongpu::float_X getV(void) const
        {
            // return velocity
            return calcVelocity(getMomentum());
        }*/
/*
        HDINLINE 
        picongpu::float_X calcGamma(const float3_X& momentum) const
        {
            //returns lorentz factor gamma = sqrt(1/(1-beta**2))
            const picongpu::float_X betaSquared = util::square<
                picongpu::float_X, 
                picongpu::float_X 
            > (calcBeta(momentum));
            return picongpu::math::sqrt( 1.0 / ( 1 - betaSquared ));
        }

        HDINLINE 
        picongpu::float_X calcBeta(const float3_X& momentum) const
        {
            // return beta = v / c
            return calcVelocity(momentum) * (1.0 / picongpu::SPEED_OF_LIGHT);
        }
*/

        HDINLINE picongpu::float_X calcBeta(const float3_X& momentum) const
        {
            // returns beta=v/c
            const picongpu::float_X gamma1 = calcGamma(momentum);
            const picongpu::float_X x = util::square<float_X, picongpu::float_X >(
                (1.0 / (mass * picongpu::SPEED_OF_LIGHT * gamma1))
            );
            return picongpu::math::sqrt((momentum * momentum).sumOfComponents() * x);
        }

        HDINLINE picongpu::float_X calcGamma(const float3_X& momentum) const
        {
            // return gamma = E/(mc^2)
            const picongpu::float_X x = util::square<float_X, picongpu::float_X > ((1.0 / (mass * picongpu::SPEED_OF_LIGHT)));
            return picongpu::math::sqrt(1.0 + (momentum * momentum).sumOfComponents() * x);
        }
/*
        HDINLINE 
        picongpu::float_X calcVelocity(const float3_X& momentum) const
        {
            //returns velocity v = p/m
           // const picongpu::float_X pOverMSquared = std::sqrt(momentum * momentum) / (mass * mass);
            const picongpu::float_X mOverPSquared = std::sqrt(mass * mass * (1.0 / (momentum * momentum)));
            //return pOverMSquared * (1.0 / (1 + pOverMSquared * (1.0 / (picongpu::SPEED_OF_LIGHT * picongpu::SPEED_OF_LIGHT))));
            return picongpu::math::sqrt(1.0/(mOverPSquared + 1.0/(picongpu::SPEED_OF_LIGHT*picongpu::SPEED_OF_LIGHT)));
        }
*/
        HDINLINE 
        picongpu::float_X calcU(const float3_X& momentum) const
        {
            //returns normalized momentum u = gama * beta
            const picongpu::float_X gamma1 = calcGamma(momentum);
            const picongpu::float_X beta1 = calcBeta(momentum);
            return gamma1 * beta1;
        }

        // Getters for Momentum in spherical coordinates
        HDINLINE
        picongpu::float_X getMomPhi(void) const
        {
            //return polar angle phi of momentum
            return calcMomPhi();
        }

        HDINLINE
        picongpu::float_X getMomTheta(void) const
        {
            //return azimuth angle psi of momentum
            return calcMomTheta();
        }

        HDINLINE
        picongpu::float_X getMomAbs(void) const
        {
            //return absolute value of momentum
            return calcMomAbs();
        }

        private:
        // Calculators for Momentum in spherical coordinates
        HDINLINE
        picongpu::float_X calcMomPhi(void) const
        {
            //return polar angle phi of momentum
            return picongpu::math::atan2(momentum.z(), momentum.x()) + picongpu::PI;
            //return picongpu::math::atan2(20.0, 1.0) + picongpu::PI;
            //return picongpu::PI / 2.0;
        }

        HDINLINE
        picongpu::float_X calcMomTheta(void) const
        {
            //return azimuth angle psi of momentum
            //return picongpu::math::acos(momentum.y() * (1.0 / 100.0));
            return picongpu::math::acos(momentum.y() * (1.0 / getMomAbs()));
            //return picongpu::math::acos(1.0 * (1.0 / 2.0));
            //return momentum.y() * 0.0;
        }

        HDINLINE
        picongpu::float_X calcMomAbs(void) const
        {
            //return absolute value of momentum
            // return picongpu::math::sqrt(momentum.x() * momentum.x()  + momentum.y() * momentum.y() + momentum.z() * momentum.z());
            return picongpu::math::sqrt((momentum * momentum).sumOfComponents());
        }
    }; // class CoolParticle
} // namespace picongpu