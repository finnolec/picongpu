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
        picongpu::float_X getCharge(void) const
        {
            return charge;
        }

        HDINLINE 
        picongpu::float_X getU(void) const
        {
            // returns normalized momentum
            return calcU(getMomentum());
        }

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

        HDINLINE
        picongpu::float_X getLocPerp(void) const
        {
            // return radial, perpendicular component of location in cylindrical coordinates
            return calcLocRho();
        }

        HDINLINE
        picongpu::float_X getLocPara(void) const
        {
            // return parallel component to z of location in cylindrical coordinates
            return location.y();
        }

        HDINLINE
        picongpu::float_X getLocPhi(void) const
        {
            //return polar angle of location in cylindrical coordinates
            return calcLocPhi();
        }

        private:
        // Calculators for Momentum in spherical coordinates
        HDINLINE
        picongpu::float_X calcMomPhi(void) const
        {
            //return polar angle phi of momentum
            return picongpu::math::atan2(momentum.z(), momentum.x()) + picongpu::PI;
        }

        HDINLINE
        picongpu::float_X calcMomTheta(void) const
        {
            //return azimuth angle psi of momentum
            //because of floating point precision x^2+y^2+z^2<y^2 for x,z<<z
            const float_X momAbs = getMomAbs();
            if(momAbs <= momentum.y())
            {
                return 0.0;
            }
            else
            {
                return picongpu::math::acos(momentum.y() * (1.0 /momAbs));
            }
        }

        HDINLINE
        picongpu::float_X calcMomAbs(void) const
        {
            //return absolute value of momentum
            // return picongpu::math::sqrt(momentum.x() * momentum.x()  + momentum.y() * momentum.y() + momentum.z() * momentum.z());
            return picongpu::math::sqrt(momentum.y() * momentum.y() + momentum.x() * momentum.x() + momentum.z() * momentum.z()); //picongpu::math::sqrt((momentum * momentum).sumOfComponents());
        }

        HDINLINE
        picongpu::float_X calcLocRho(void) const
        {
            // return radial component of location in cylindrical coordinates
            return picongpu::math::sqrt(location.x() * location.x() + location.z() * location.z());
        }

        HDINLINE
        picongpu::float_X calcLocPhi(void) const
        {
            // return radial angle phi of location in cylindrical coordinates
            return picongpu::math::atan2(location.z(), location.x()) + picongpu::PI;
        }

    }; // class CoolParticle
} // namespace picongpu