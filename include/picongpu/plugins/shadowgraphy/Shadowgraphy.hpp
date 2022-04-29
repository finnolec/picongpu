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
#include <pmacc/mpi/MPIReduce.hpp>
#include <pmacc/mpi/reduceMethods/Reduce.hpp>

#include <sstream>
#include <iostream>
#include <string>

#include <stdio.h>

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

#include "picongpu/plugins/shadowgraphy/ShadowgraphyHelper.hpp"

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

                std::unique_ptr<container::DeviceBuffer<float3_64, 2>> dBuffer_SI1;
                std::unique_ptr<container::DeviceBuffer<float3_64, 2>> dBuffer_SI2;

                bool isIntegrating;
                int startTime;

                bool isMaster = false;

                shadowgraphy::Helper* helper = nullptr;
                pmacc::mpi::MPIReduce reduce;
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

                        // Time integration from param files
                        std::cout<<"hello world what what what "<< SI::DELTA_T_SI << std::endl;
                        std::cout<<std::stoi(this->notifyPeriod) << std::endl;
                        int startTime = std::stoi(this->notifyPeriod);
                        int endTime = std::stoi(this->notifyPeriod) + params::t_n;

                        std::cout<<endTime<<std::endl;
                        std::string internalNotifyPeriod = std::to_string(startTime) + ":" + std::to_string(endTime) + ":" + std::to_string(params::t_res);

                        std::cout<<internalNotifyPeriod<<std::endl;
                        
                        Environment<>::get().PluginConnector().setNotificationPeriod(this, internalNotifyPeriod);
                        namespace vec = ::pmacc::math;
                        typedef SuperCellSize BlockDim;

                        vec::Size_t<simDim> size = vec::Size_t<simDim>(this->cellDescription->getGridSuperCells())
                                * precisionCast<size_t>(BlockDim::toRT())
                            - precisionCast<size_t>(2 * BlockDim::toRT());
                        this->dBuffer_SI1 = std::make_unique<container::DeviceBuffer<float3_64, simDim - 1>>(
                            size.shrink<simDim - 1>((this->plane + 1) % simDim));
                        this->dBuffer_SI2 = std::make_unique<container::DeviceBuffer<float3_64, simDim - 1>>(
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
                    //std::cout << "Shadowgraphy notify period is: " << currentStep << std::endl;
                    //std::cout << "231 " << picongpu::SI::DELTA_T_SI << std::endl;
                    //printf("231 %e\n", float(picongpu::SI::DELTA_T_SI));

                    if(sliceIsOK)
                    {
                        
                        isMaster = reduce.hasResult(pmacc::mpi::reduceMethods::Reduce());

                        // First time the plugin is called:
                        if(isIntegrating == false)
                        {
                            if (isMaster)
                            {
                                // Get grid size
                                namespace vec = pmacc::math;
                                typedef SuperCellSize BlockDim;
                                DataConnector& dc = Environment<>::get().DataConnector();
                                auto field = dc.get<FieldE>(FieldE::getName(), true)
                                                            ->getGridBuffer()
                                                            .getDeviceBuffer()
                                                            .cartBuffer()
                                                            .view(BlockDim::toRT(), -BlockDim::toRT());

                                pmacc::GridController<simDim>& con = pmacc::Environment<simDim>::get().GridController();
                                vec::Size_t<simDim> gpuDim = (vec::Size_t<simDim>) con.getGpuNodes();
                                vec::Size_t<simDim> globalGridSize = gpuDim * field.size();

                                helper = new Helper(globalGridSize);
                            }


                            // Create Integrator object %TODO
                            startTime = currentStep;
                            isIntegrating = true;
                        }

                        int localStep = (currentStep - startTime) / params::t_res;

                        std::cout << "localStep: " << localStep << std::endl;

                        if(localStep != int(params::t_n / params::t_res))
                        {
                            namespace vec = ::pmacc::math;
                            typedef SuperCellSize BlockDim;
                            DataConnector& dc = Environment<>::get().DataConnector();
                            auto field_coreBorderE = dc.get<FieldE>(FieldE::getName(), true)
                                                        ->getGridBuffer()
                                                        .getDeviceBuffer()
                                                        .cartBuffer()
                                                        .view(BlockDim::toRT(), -BlockDim::toRT());

                            storeSlice<FieldE>(field_coreBorderE, this->plane, this->slicePoint, localStep);

                            auto field_coreBorderB = dc.get<FieldB>(FieldB::getName(), true)
                                                        ->getGridBuffer()
                                                        .getDeviceBuffer()
                                                        .cartBuffer()
                                                        .view(BlockDim::toRT(), -BlockDim::toRT());

                            storeSlice<FieldB>(field_coreBorderB, this->plane, this->slicePoint, localStep);

                            helper->calculate_energy_flux(localStep, true);
                            helper->calculate_energy_flux(localStep, false);
                        }
                        else
                        {
                            std::ostringstream filename;
                            filename << this->fileName << "_" << startTime << ":" << currentStep << ".dat";

                            //data = helper->get_shadowgram();
                            writeFile(helper->get_shadowgram(), filename.str());

                            std::cout << "destructor called" << std::endl;
                            delete(helper);
                            isIntegrating = false;
                        }
                    }
                }

                template<typename Field, typename TField>
                void storeSlice(const TField& field, int nAxis, float slicePoint, int localStep)
                {
                    namespace vec = pmacc::math;

                    pmacc::GridController<simDim>& con = pmacc::Environment<simDim>::get().GridController();
                    vec::Size_t<simDim> gpuDim = (vec::Size_t<simDim>) con.getGpuNodes();
                    vec::Size_t<simDim> globalGridSize = gpuDim * field.size();

                    // FIRST SLICE OF FIELD FOR YEE OFFSET
                    int globalPlane1 = globalGridSize[nAxis] * slicePoint;
                    int localPlane1 = globalPlane1 % field.size()[nAxis];
                    int gpuPlane1 = globalPlane1 / field.size()[nAxis];

                    vec::Int<simDim> nVector1(vec::Int<simDim>::create(0));
                    nVector1[nAxis] = 1;

                    zone::SphericZone<simDim> gpuGatheringZone1(gpuDim, nVector1 * gpuPlane1);
                    gpuGatheringZone1.size[nAxis] = 1;

                    algorithm::mpi::Gather<simDim> gather(gpuGatheringZone1);

                    if(!gather.participate())
                        return;

                    vec::UInt32<3> twistedAxesVec1((nAxis + 1) % 3, (nAxis + 2) % 3, nAxis);

                    // convert data to higher precision and to SI units
                    ShadowgraphyHelper::ConversionFunctor<Field> cf1;
                    algorithm::kernel::RT::Foreach()(
                        dBuffer_SI1->zone(),
                        dBuffer_SI1->origin(),
                        cursor::tools::slice(field.originCustomAxes(twistedAxesVec1)(0, 0, localPlane1)),
                        cf1);
            

                    // copy selected plane from device to host
                    container::HostBuffer<float3_64, simDim - 1> hBuffer1(dBuffer_SI1->size());
                    hBuffer1 = *dBuffer_SI1;

                    // collect data from all nodes/GPUs
                    vec::Size_t<simDim> globalDomainSize1 = Environment<simDim>::get().SubGrid().getGlobalDomain().size;
                    vec::Size_t<simDim - 1> globalSliceSize1 = globalDomainSize1.shrink<simDim - 1>((nAxis + 1) % simDim);
                    container::HostBuffer<float3_64, simDim - 1> globalBuffer1(globalSliceSize1);
                    gather(globalBuffer1, hBuffer1, nAxis);
                    if(!gather.root())
                        return;

                    
                    // SECOND SLICE OF FIELD FOR YEE OFFSET
                    int globalPlane2 = globalGridSize[nAxis] * slicePoint + 1;
                    int localPlane2 = globalPlane2 % field.size()[nAxis];
                    int gpuPlane2 = globalPlane2 / field.size()[nAxis];

                    vec::Int<simDim> nVector2(vec::Int<simDim>::create(0));
                    nVector2[nAxis] = 1;

                    zone::SphericZone<simDim> gpuGatheringZone2(gpuDim, nVector2 * gpuPlane2);
                    gpuGatheringZone2.size[nAxis] = 1;

                    //algorithm::mpi::Gather<simDim> gather(gpuGatheringZone2);

                    if(!gather.participate())
                        return;

                    vec::UInt32<3> twistedAxesVec2((nAxis + 1) % 3, (nAxis + 2) % 3, nAxis);

                    // convert data to higher precision and to SI units 
                    ShadowgraphyHelper::ConversionFunctor<Field> cf2;
                    algorithm::kernel::RT::Foreach()(
                        dBuffer_SI2->zone(),
                        dBuffer_SI2->origin(),
                        cursor::tools::slice(field.originCustomAxes(twistedAxesVec2)(0, 0, localPlane2)),
                        cf2);
            

                    // copy selected plane from device to host 
                    container::HostBuffer<float3_64, simDim - 1> hBuffer2(dBuffer_SI2->size());
                    hBuffer2 = *dBuffer_SI2;

                    /// collect data from all nodes/GPUs 
                    vec::Size_t<simDim> globalDomainSize2 = Environment<simDim>::get().SubGrid().getGlobalDomain().size;
                    vec::Size_t<simDim - 1> globalSliceSize2 = globalDomainSize2.shrink<simDim - 1>((nAxis + 1) % simDim);
                    container::HostBuffer<float3_64, simDim - 1> globalBuffer2(globalSliceSize2);
                    gather(globalBuffer2, hBuffer2, nAxis);
                    if(!gather.root())
                        return;
                    

                    if(isMaster)
                    {
                        helper->store_field<Field>(localStep, &globalBuffer1, &globalBuffer2);
                    }
                    //std::ofstream file(filename.c_str());
                    //file << globalBuffer;

                }

                void writeFile(std::vector< std::vector< float_64 > > values, std::string name)
                {
                    std::ofstream outFile;
                    outFile.open(name.c_str(), std::ofstream::out | std::ostream::trunc);

                    if(!outFile)
                    {
                        std::cerr << "Can't open file [" << name << "] for output, disable plugin output. "
                                  << std::endl;
                        isMaster = false; // no Master anymore -> no process is able to write
                    }
                    else
                    {
                        for( unsigned int i = 0; i < helper->get_n_x(); ++i ) // over all x
                        {
                            for(unsigned int j = 0;  j < helper->get_n_y(); ++j) // over all y
                            {
                                outFile << values[i][j] << "\t";
                            } // for loop over all y

                            outFile << std::endl;
                        } // for loop over all x

                        outFile.flush();
                        outFile << std::endl; // now all data are written to file

                        if(outFile.fail())
                            std::cerr << "Error on flushing file [" << name << "]. " << std::endl;

                        outFile.close();
                    }
                }
            };
        }
    }
}