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

#include "picongpu/fields/FieldB.hpp"
#include "picongpu/fields/FieldE.hpp"

#include <pmacc/cuSTL/algorithm/host/Foreach.hpp>
#include <pmacc/cuSTL/algorithm/kernel/run-time/Foreach.hpp>
#include <pmacc/cuSTL/algorithm/mpi/Gather.hpp>
#include <pmacc/cuSTL/container/DeviceBuffer.hpp>
#include <pmacc/cuSTL/container/HostBuffer.hpp>
#include <pmacc/cuSTL/cursor/tools/slice.hpp>
#include <pmacc/dataManagement/DataConnector.hpp>
#include <pmacc/math/Vector.hpp>
#include <pmacc/math/vector/Float.hpp>
#include <pmacc/math/vector/Int.hpp>
#include <pmacc/math/vector/Size_t.hpp>

#include <sstream>
#include <iostream>
#include <string>

#include "picongpu/fields/FieldB.hpp"
#include "picongpu/fields/FieldE.hpp"

#include <pmacc/cuSTL/algorithm/host/Foreach.hpp>
#include <pmacc/cuSTL/algorithm/kernel/Foreach.hpp>
#include <pmacc/cuSTL/algorithm/mpi/Gather.hpp>
#include <pmacc/cuSTL/container/DeviceBuffer.hpp>
#include <pmacc/cuSTL/container/HostBuffer.hpp>
#include <pmacc/cuSTL/cursor/tools/slice.hpp>
#include <pmacc/dataManagement/DataConnector.hpp>
#include <pmacc/math/Vector.hpp>
#include <pmacc/math/vector/Float.hpp>
#include <pmacc/math/vector/Int.hpp>
#include <pmacc/math/vector/Size_t.hpp>

#include <sstream>

namespace picongpu
{
    using namespace pmacc;
    namespace po = boost::program_options;

    namespace plugins
    {
        namespace shadowgraphy
        {

            namespace ShadowgraphyHelper
            {
                template<class Field>
                class ConversionFunctor
                {
                public:
                    /* convert field data to higher precision and convert to SI units on GPUs */
                    template<typename T_Acc>
                    DINLINE void operator()(T_Acc const& acc, float3_64& target, const typename Field::ValueType fieldData)
                        const
                    {
                        target = precisionCast<float_64>(fieldData) * float_64((Field::getUnit())[0]);
                    }
                };
            } // end namespace ShadowgraphyHelper

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
                std::string fileName;
                float_X slicePoint;

                std::unique_ptr<container::DeviceBuffer<float3_64, 2>> dBuffer_SI;

                bool isIntegrating;
                int starttime;
    /*
            std::string name;
            std::string prefix;
            std::vector<std::string> notifyPeriod;
            std::vector<std::string> fileName;
            std::vector<int> plane;
            std::vector<float_X> slicePoint;
            MappingDesc* cellDescription;
            std::vector<SliceFieldPrinter<Field>> childs;
    */

            public:
                Shadowgraphy()
                    : pluginName("TransitionRadiation: calculate transition radiation of species")
                    , isIntegrating(false)
                {
                    /* register our plugin during creation */
                    Environment<>::get().PluginConnector().registerPlugin(this);
                    pluginPrefix = "shadowgraphy";
                }

                std::string pluginGetName() const override
                {
                    return "Shadowgraphy";
                }


                void pluginRegisterHelp(po::options_description& desc) override
                {
                    /* register command line parameters for your plugin */
                    //desc.add_options()
                    //("shadowgraphy.period", po::value<uint32_t > (&notifyPeriod)->default_value(0),
                    //"Enable Shadowgraphy [for each n-th step]")   
                    desc.add_options()(
                        (this->pluginPrefix + ".period").c_str(),
                        po::value<std::string>(&this->notifyPeriod)->multitoken(),
                        "notify period");
                    desc.add_options()(
                        (this->pluginPrefix + ".fileName").c_str(),
                        po::value<std::string>(&this->fileName)->multitoken(),
                        "file name to store slices in");
                    desc.add_options()(
                        (this->pluginPrefix + ".plane").c_str(),
                        po::value<int>(&this->plane)->multitoken(),
                        "specifies the axis which stands on the cutting plane (0,1,2)");
                    desc.add_options()(
                        (this->pluginPrefix + ".slicePoint").c_str(),
                        po::value<float_X>(&this->slicePoint)->multitoken(),
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
                        //this->dBBuffer_SI = std::make_unique<container::DeviceBuffer<float3_64, simDim - 1>>(
                        //    size.shrink<simDim - 1>((this->plane + 1) % simDim));
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


                void notify(uint32_t currentStep) override
                {
                    /* notification callback for simulation step currentStep
                    * called every notifyPeriod steps */
                    std::cout << "Shadowgraphy notify period is: " << currentStep;


                    if(sliceIsOK)
                    {
                        // First time the plugin is called:
                        if(isIntegrating == false)
                        {
                            // Create Integrator object %TODO
                            startTime = currentStep;
                            isIntegrating = true;
                        }

                        int localstep = currenStep - startTime;

                        if(localstep != params::t_n)
                        {
                            namespace vec = ::pmacc::math;
                            typedef SuperCellSize BlockDim;
                            DataConnector& dc = Environment<>::get().DataConnector();
                            auto field_coreBorder = dc.get<FieldE>(FieldE::getName(), true)
                                                        ->getGridBuffer()
                                                        .getDeviceBuffer()
                                                        .cartBuffer()
                                                        .view(BlockDim::toRT(), -BlockDim::toRT());

                            std::ostringstream filenameE;
                            filenameE << this->fileName << "_E_" << currentStep << ".dat";
                            storeSlice<FieldE>(field_coreBorder, this->plane, this->slicePoint, filenameE.str());

                            field_coreBorder = dc.get<FieldB>(FieldB::getName(), true)
                                                        ->getGridBuffer()
                                                        .getDeviceBuffer()
                                                        .cartBuffer()
                                                        .view(BlockDim::toRT(), -BlockDim::toRT());

                            std::ostringstream filenameB;
                            filenameB << this->fileName << "_B_" << currentStep << ".dat";
                            storeSlice<FieldB>(field_coreBorder, this->plane, this->slicePoint, filenameB.str());
                            
                        }
                        else
                        {
                            //delete(Integrator) %TODO
                            isIntegrating = false;
                        }
                        
                    }
                }

                template<typename Field, typename TField>
                void storeSlice(const TField& field, int nAxis, float slicePoint, std::string filename)
                {
                    namespace vec = pmacc::math;

                    pmacc::GridController<simDim>& con = pmacc::Environment<simDim>::get().GridController();
                    vec::Size_t<simDim> gpuDim = (vec::Size_t<simDim>) con.getGpuNodes();
                    vec::Size_t<simDim> globalGridSize = gpuDim * field.size();
                    int globalPlane = globalGridSize[nAxis] * slicePoint;
                    int localPlane = globalPlane % field.size()[nAxis];
                    int gpuPlane = globalPlane / field.size()[nAxis];

                    vec::Int<simDim> nVector(vec::Int<simDim>::create(0));
                    nVector[nAxis] = 1;

                    zone::SphericZone<simDim> gpuGatheringZone(gpuDim, nVector * gpuPlane);
                    gpuGatheringZone.size[nAxis] = 1;

                    algorithm::mpi::Gather<simDim> gather(gpuGatheringZone);

                    if(!gather.participate())
                        return;

                    vec::UInt32<3> twistedAxesVec((nAxis + 1) % 3, (nAxis + 2) % 3, nAxis);

                    /* convert data to higher precision and to SI units */
                    ShadowgraphyHelper::ConversionFunctor<Field> cf;
                    algorithm::kernel::RT::Foreach()(
                        dBuffer_SI->zone(),
                        dBuffer_SI->origin(),
                        cursor::tools::slice(field.originCustomAxes(twistedAxesVec)(0, 0, localPlane)),
                        cf);
            

                    /* copy selected plane from device to host */
                    container::HostBuffer<float3_64, simDim - 1> hBuffer(dBuffer_SI->size());
                    hBuffer = *dBuffer_SI;

                    /* collect data from all nodes/GPUs */
                    vec::Size_t<simDim> globalDomainSize = Environment<simDim>::get().SubGrid().getGlobalDomain().size;
                    vec::Size_t<simDim - 1> globalSliceSize = globalDomainSize.shrink<simDim - 1>((nAxis + 1) % simDim);
                    container::HostBuffer<float3_64, simDim - 1> globalBuffer(globalSliceSize);
                    gather(globalBuffer, hBuffer, nAxis);
                    if(!gather.root())
                        return;

                    std::ofstream file(filename.c_str());
                    file << globalBuffer;
                }
            };
        }
    }
}