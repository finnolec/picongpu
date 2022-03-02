/* Copyright 2013-2022 Axel Huebl, Heiko Burau, Rene Widera, Richard Pausch,
 *                     Klaus Steiniger, Felix Schmitt, Benjamin Worpitz
 *                     Finn-Ole Carstens
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

#include "picongpu/plugins/ILightweightPlugin.hpp"

#include <iostream>

namespace picongpu
{
    namespace plugins
    {

        class Shadowgraphy : public ILightweightPlugin
        {
        private:
            // technical variables for PIConGPU plugins
            std::string pluginName;
            std::string pluginPrefix;
            
            MappingDesc* cellDescription = nullptr;
            std::string notifyPeriod;

            bool sliceIsOK;
            int plane;


        public:
            Shadowgraphy()
                : pluginName("TransitionRadiation: calculate transition radiation of species")
            {
                /* register our plugin during creation */
                Environment<>::get().PluginConnector().registerPlugin(this);
                pluginPrefix = "shadowgraphy";
            }

            std::string pluginGetName() const override
            {
                return "Shadowgraphy";
            }

            void notify(uint32_t currentStep) override
            {
                /* notification callback for simulation step currentStep
                * called every notifyPeriod steps */
                std::cout << "Shadowgraphy notify period is: " << currentStep;

                resetBuffers();

                calculateTransitionRadiation(currentStep);

                collectDataGPUToMaster();
            }

            void pluginRegisterHelp(po::options_description& desc) override
            {
                /* register command line parameters for your plugin */
                //desc.add_options()
                //("shadowgraphy.period", po::value<uint32_t > (&notifyPeriod)->default_value(0),
                //"Enable Shadowgraphy [for each n-th step]");
                desc.add_options()(
                    (pluginPrefix + ".period").c_str(),
                    po::value<std::string>(&notifyPeriod),
                    "enable plugin [for each n-th step]");
                desc.add_options()(
                    (this->prefix + ".plane").c_str(),
                    po::value<std::vector<int>>(&this->plane)->multitoken(),
                    "specifies the axis which stands on the cutting plane (0,1,2)");
                desc.add_options()(
                    (pluginPrefix + ".slicePoint").c_str(),
                    po::value<std::vector<float_X>>(&this->slicePoint)->multitoken(),
                    "slice point 0.0 <= x <= 1.0");
            }

            void pluginLoad() override
            {
                /* called when plugin is loaded, command line flags are available here
                * set notification period for our plugin at the PluginConnector */
                if(float_X(0.0) <= slicePoint && slicePoint <= float_X(1.0))
                {
                    /* in case the slice point is inside of [0.0,1.0] */
                    sliceIsOK = true;
                    Environment<>::get().PluginConnector().setNotificationPeriod(this, this->notifyPeriod);
                    namespace vec = ::pmacc::math;
                    typedef SuperCellSize BlockDim;

                    vec::Size_t<simDim> size = vec::Size_t<simDim>(this->cellDescription->getGridSuperCells())
                            * precisionCast<size_t>(BlockDim::toRT())
                        - precisionCast<size_t>(2 * BlockDim::toRT());
                    this->dBuffer_SI = std::make_unique<container::DeviceBuffer<float3_64, simDim - 1>>(
                        size.shrink<simDim - 1>((this->plane + 1) % simDim));
                }
                else
                {
                    /* in case the slice point is outside of [0.0,1.0] */
                    sliceIsOK = false;
                    std::cerr << "In the SliceFieldPrinter plugin a slice point"
                            << " (slice_point=" << slicePoint << ") is outside of [0.0, 1.0]. " << std::endl
                            << "The request will be ignored. " << std::endl;
                }
            }

            void pluginUnload() override
            {
                /* called when plugin is unloaded, cleanup here */
            }

            /** Implementation of base class function. Sets mapping description.
                *
                * @param cellDescription
                */
            void setMappingDescription(MappingDesc* cellDescription) override
            {
                this->cellDescription = cellDescription;
            }

        private:
            //! Moves transition radiation data from GPUs to CPUs.
            void copyRadiationDeviceToHost()
            {
                incTransRad->deviceToHost();
                __getTransactionEvent().waitForFinished();
                cohTransRadPara->deviceToHost();
                __getTransactionEvent().waitForFinished();
                cohTransRadPerp->deviceToHost();
                __getTransactionEvent().waitForFinished();
                numParticles->deviceToHost();
                __getTransactionEvent().waitForFinished();
            }

            /** Combine transition radiation data from each CPU and store result on master.
            *
            * @remark copyRadiationDeviceToHost( ) should be called before.
            */
            void collectRadiationOnMaster()
            {
                reduce(
                    pmacc::math::operation::Add(),
                    tmpITR.data(),
                    incTransRad->getHostBuffer().getBasePointer(),
                    elementsTransitionRadiation(),
                    mpi::reduceMethods::Reduce());
                reduce(
                    pmacc::math::operation::Add(),
                    tmpCTRpara.data(),
                    cohTransRadPara->getHostBuffer().getBasePointer(),
                    elementsTransitionRadiation(),
                    mpi::reduceMethods::Reduce());
                reduce(
                    pmacc::math::operation::Add(),
                    tmpCTRperp.data(),
                    cohTransRadPerp->getHostBuffer().getBasePointer(),
                    elementsTransitionRadiation(),
                    mpi::reduceMethods::Reduce());
                reduce(
                    pmacc::math::operation::Add(),
                    tmpNum.data(),
                    numParticles->getHostBuffer().getBasePointer(),
                    elementsTransitionRadiation(),
                    mpi::reduceMethods::Reduce());
            }

            //! perform all operations to get data from GPU to master
            void collectDataGPUToMaster()
            {
                // collect data GPU -> CPU -> Master
                copyRadiationDeviceToHost();
                collectRadiationOnMaster();
                sumTransitionRadiation(
                    theTransRad.data(),
                    tmpITR.data(),
                    tmpCTRpara.data(),
                    tmpCTRperp.data(),
                    tmpNum.data());
            }
        };
    }
}