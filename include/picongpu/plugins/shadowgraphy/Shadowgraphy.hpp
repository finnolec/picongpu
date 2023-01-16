/* Copyright 2023 Finn-Ole Carstens, Rene Widera
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

#include "picongpu/fields/FieldB.hpp"
#include "picongpu/fields/FieldE.hpp"
#include "picongpu/plugins/ILightweightPlugin.hpp"
#include "picongpu/plugins/common/openPMDAttributes.hpp"
#include "picongpu/plugins/common/openPMDDefaultExtension.hpp"
#include "picongpu/plugins/common/openPMDVersion.def"
#include "picongpu/plugins/common/openPMDWriteMeta.hpp"
#include "picongpu/plugins/shadowgraphy/GatherSlice.hpp"
#include "picongpu/plugins/shadowgraphy/ShadowgraphyHelper.hpp"

#include <pmacc/dataManagement/DataConnector.hpp>
#include <pmacc/math/Vector.hpp>

#include <iostream>
#include <sstream>
#include <string>

#include <mpi.h>
#include <openPMD/openPMD.hpp>
#include <stdio.h>


namespace picongpu
{
    using namespace pmacc;
    namespace po = boost::program_options;

    namespace plugins
    {
        namespace shadowgraphy
        {
            class Shadowgraphy : public ILightweightPlugin
            {
            private:
                // technical variables for PIConGPU plugins
                std::string pluginName;
                std::string pluginPrefix;
                std::string filenameExtension = openPMD::getDefaultExtension().c_str();

                MappingDesc* m_cellDescription = nullptr;
                std::string notifyPeriod;

                bool sliceIsOK = false;

                // do not change the plane, code is only supporting a plane in z direction
                int plane = 2;
                std::string fileName;
                float_X slicePoint = 0.5;

                bool isIntegrating = false;
                int startTime = 0;

                ///@todo: we must verify that this command line param is set!
                float_X focuspos = 0.0_X;
                int duration = 0;

                int localPlaneIdx = -1;

                std::unique_ptr<shadowgraphy::Helper> helper;
                std::unique_ptr<shadowgraphy::GatherSlice> gather;

                bool fourierOutputEnabled = false;
                bool intermediateOutputEnabled = false;

            public:
                Shadowgraphy()
                    : pluginName("Shadowgraphy: calculate the energy density of a laser by integrating"
                                 "the Poynting vectors in a spatially fixed slice over a given time interval")

                {
                    static_assert(simDim == DIM3, "Shadowgraphy-plugin requires 3D simulations.");

                    /* register our plugin during creation */
                    Environment<>::get().PluginConnector().registerPlugin(this);
                    pluginPrefix = "shadowgraphy";
                }

                //! Implementation of base class function.
                std::string pluginGetName() const override
                {
                    return "Shadowgraphy";
                }


                //! Implementation of base class function.
                void pluginRegisterHelp(po::options_description& desc) override
                {
#if(PIC_ENABLE_FFTW3 == 1)
                    desc.add_options()(
                        (this->pluginPrefix + ".period").c_str(),
                        po::value<std::string>(&this->notifyPeriod)->multitoken(),
                        "notify period");
                    desc.add_options()(
                        (this->pluginPrefix + ".fileName").c_str(),
                        po::value<std::string>(&this->fileName)->multitoken(),
                        "file name to store slices in");
                    desc.add_options()(
                        (this->pluginPrefix + ".ext").c_str(),
                        po::value<std::string>(&this->filenameExtension)->multitoken(),
                        "openPMD filename extension");
                    desc.add_options()(
                        (this->pluginPrefix + ".slicePoint").c_str(),
                        po::value<float_X>(&this->slicePoint)->multitoken(),
                        "slice point 0.0 <= x < 1.0");
                    desc.add_options()(
                        (this->pluginPrefix + ".focuspos").c_str(),
                        po::value<float_X>(&this->focuspos)->multitoken(),
                        "focus position relative to slice point in microns");
                    desc.add_options()(
                        (this->pluginPrefix + ".duration").c_str(),
                        po::value<int>(&this->duration)->multitoken(),
                        "nt");
                    desc.add_options()(
                        (this->pluginPrefix + ".fourieroutput").c_str(),
                        po::value<bool>(&fourierOutputEnabled)->zero_tokens(),
                        "optional output: E and B fields in (kx, ky, omega) Fourier space");
                    desc.add_options()(
                        (this->pluginPrefix + ".intermediateoutput").c_str(),
                        po::value<bool>(&intermediateOutputEnabled)->zero_tokens(),
                        "optional output: E and B fields in (x, y, omega) Fourier space");
#else
                    desc.add_options()(
                        (this->pluginPrefix).c_str(),
                        "plugin disabled [compiled without dependency FFTW]");
#endif
                }
                plugins::multi::Option<std::string> extension
                    = {"ext", "openPMD filename extension", openPMD::getDefaultExtension().c_str()};

                //! Implementation of base class function.
                void pluginLoad() override
                {
                    /* called when plugin is loaded, command line flags are available here
                     * set notification period for our plugin at the PluginConnector */
                    if(0 != notifyPeriod.size())
                    {
                        if(float_X(0.0) <= slicePoint && slicePoint < float_X(1.0))
                        {
                            /* in case the slice point is inside of [0.0,1.0) */
                            sliceIsOK = true;

                            /* The plugin integrates the Poynting vectors over time and must thus be called every
                             * tRes-th time-step of the simulation until the integration is done */
                            int startTime = std::stoi(this->notifyPeriod);
                            int endTime = std::stoi(this->notifyPeriod) + this->duration;

                            std::string internalNotifyPeriod = std::to_string(startTime) + ":"
                                + std::to_string(endTime) + ":" + std::to_string(params::tRes);

                            Environment<>::get().PluginConnector().setNotificationPeriod(this, internalNotifyPeriod);

                            const SubGrid<simDim>& subGrid = Environment<simDim>::get().SubGrid();
                            auto globalDomain = subGrid.getGlobalDomain();
                            auto globalPlaneExtent = globalDomain.size[plane];
                            auto localDomain = subGrid.getLocalDomain();

                            auto globalPlaneIdx = globalPlaneExtent * slicePoint;

                            auto isPlaneInLocalDomain = globalPlaneIdx >= localDomain.offset[plane]
                                && globalPlaneIdx < localDomain.offset[plane] + localDomain.size[plane];
                            if(isPlaneInLocalDomain)
                                localPlaneIdx = globalPlaneIdx - localDomain.offset[plane];

                            gather = std::make_unique<shadowgraphy::GatherSlice>();
                            gather->participate(isPlaneInLocalDomain);
                        }
                        else
                        {
                            /* in case the slice point is outside of [0.0,1.0) */
                            sliceIsOK = false;
                            std::cerr << "In the Shadowgraphy plugin the slice point"
                                      << " (slicePoint=" << slicePoint << ") is outside of [0.0, 1.0). " << std::endl
                                      << "The request will be ignored. " << std::endl;
                        }
                    }
                    else
                    {
                        sliceIsOK = false;
                    }
                }

                //! Implementation of base class function.
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
                    this->m_cellDescription = cellDescription;
                }

                //! Implementation of base class function.
                void notify(uint32_t currentStep) override
                {
                    // skip notify, slice is not intersecting the local domain
                    if(!gather->isParticipating())
                        return;
                    /* notification callback for simulation step currentStep
                     * called every notifyPeriod steps */
                    if(sliceIsOK)
                    {
                        // First time the plugin is called:
                        if(isIntegrating == false)
                        {
                            startTime = currentStep;

                            if(gather->isMaster() && helper == nullptr)
                            {
                                std::cout << "master init " << currentStep << std::endl;
                                /// @todo shared pointer
                                helper = std::make_unique<Helper>(
                                    currentStep,
                                    this->slicePoint,
                                    this->focuspos * 1e-6,
                                    this->duration,
                                    this->fourierOutputEnabled,
                                    this->intermediateOutputEnabled);
                            }
                            // Create Integrator object %TODO
                            isIntegrating = true;
                        }

                        // convert currentStep (simulation time-step) into localStep for time domain DFT
                        int localStep = (currentStep - startTime) / params::tRes;

                        std::cout << "try calculate " << localStep << std::endl;

                        bool const dumpFinalData = localStep == int(this->duration / params::tRes);
                        if(!dumpFinalData)
                        {
                            std::cout << "prepare E field" << std::endl;
                            auto sliceBufferE
                                = getGlobalSlice<shadowgraphy::Helper::FieldType::E>(FieldE::getName(), localPlaneIdx);
                            if(gather->isMaster())
                            {
                                std::cout << " finish preparing global slice" << std::endl;
                                helper->storeField<shadowgraphy::Helper::FieldType::E>(
                                    localStep,
                                    currentStep,
                                    sliceBufferE);
                            }

                            std::cout << "prepare B field" << std::endl;
                            auto sliceBufferB
                                = getGlobalSlice<shadowgraphy::Helper::FieldType::B>(FieldB::getName(), localPlaneIdx);
                            if(gather->isMaster())
                            {
                                std::cout << " finish preparing global slice" << std::endl;
                                helper->storeField<shadowgraphy::Helper::FieldType::B>(
                                    localStep,
                                    currentStep,
                                    sliceBufferB);
                            }

                            if(gather->isMaster())
                            {
                                helper->calculate_dft(localStep);
                            }
                        }
                        else
                        {
                            if(gather->isMaster())
                            {
                                std::cout << "dump " << currentStep << std::endl;
                                helper->propagateFieldsAndCalculateShadowgram();

                                std::ostringstream filename;
                                filename << this->fileName << "_" << startTime << ":" << currentStep << ".dat";

                                writeFile(helper->getShadowgram(), filename.str());

                                writeToOpenPMDFile(currentStep);

                                // delete helper and free all memory
                                helper.reset(nullptr);
                            }
                            isIntegrating = false;
                        }
                    }
                }

            private:
                /** Create and store the global slice out of local data.
                 *
                 * Create the slice of the local field. The field values will be interpolated to the origin of the
                 * cell. Gather the local field data into a single global field on the gather master.
                 *
                 * @tparam T_fieldType
                 * @param fieldName
                 * @param localPlaneIdx
                 * @return Buffer with gathered global slice. (only gather master buffer contains data)
                 */
                template<typename shadowgraphy::Helper::FieldType T_fieldType>
                auto getGlobalSlice(std::string fieldName, int localPlaneIdx) const
                    -> std::shared_ptr<HostBufferIntern<float2_X, DIM2>>
                {
                    DataConnector& dc = Environment<>::get().DataConnector();
                    auto inputFieldBuffer = dc.get<FieldE>(FieldE::getName(), false);

                    const SubGrid<simDim>& subGrid = Environment<simDim>::get().SubGrid();
                    auto globalDomain = subGrid.getGlobalDomain();
                    auto globalPlaneExtent = globalDomain.size[plane];

                    auto localDomainOffset = subGrid.getLocalDomain().offset.shrink<DIM2>(0);
                    auto globalDomainSliceSize = subGrid.getGlobalDomain().size.shrink<DIM2>(0);

                    auto fieldSlice = createSlice<T_fieldType>(inputFieldBuffer, localPlaneIdx);
                    return gather->gatherSlice(fieldSlice, globalDomainSliceSize, localDomainOffset);
                }

                template<typename shadowgraphy::Helper::FieldType T_fieldType, typename T_FieldBuffer>
                auto createSlice(std::shared_ptr<T_FieldBuffer> inputFieldBuffer, int sliceCellZ) const
                {
                    auto bufferGridLayout = inputFieldBuffer->getGridLayout();
                    DataSpace<DIM2> localSliceSize
                        = bufferGridLayout.getDataSpaceWithoutGuarding().template shrink<DIM2>(0);

                    // skip guard cells
                    auto inputFieldBox = inputFieldBuffer->getHostDataBox().shift(bufferGridLayout.getGuard());

                    auto sliceBuffer = std::make_shared<HostBufferIntern<float2_X, DIM2>>(localSliceSize);
                    auto sliceBox = sliceBuffer->getDataBox();

                    std::cout << " start loading slice" << std::endl;
                    for(int y = 0; y < localSliceSize.y(); ++y)
                        for(int x = 0; x < localSliceSize.x(); ++x)
                        {
                            DataSpace<DIM2> idx(x, y);
                            DataSpace<DIM3> srcIdx(idx.x(), idx.y(), sliceCellZ);
                            sliceBox(idx) = helper->cross<T_fieldType>(inputFieldBox.shift(srcIdx));
                        }
                    std::cout << " end loading slice" << std::endl;

                    return sliceBuffer;
                }

                void writeToOpenPMDFile(uint32_t currentStep)
                {
                    std::stringstream filename;
                    filename << pluginPrefix << "_%T." << filenameExtension;
                    ::openPMD::Series series(filename.str(), ::openPMD::Access::CREATE);

                    ::openPMD::Extent extent
                        = {static_cast<unsigned long int>(helper->getSizeY()),
                           static_cast<unsigned long int>(helper->getSizeX())};
                    ::openPMD::Offset offset = {0, 0};
                    ::openPMD::Datatype datatype = ::openPMD::determineDatatype<float_64>();
                    ::openPMD::Dataset dataset{datatype, extent};

                    auto mesh = series.iterations[currentStep].meshes["shadowgram"];
                    mesh.setAxisLabels(std::vector<std::string>{"x", "y"});
                    mesh.setDataOrder(::openPMD::Mesh::DataOrder::F);
                    mesh.setGridUnitSI(1.0);
                    mesh.setGridSpacing(std::vector<double>{1.0, 1.0});
                    mesh.setGeometry(::openPMD::Mesh::Geometry::cartesian); // set be default

                    auto shadowgram = mesh[::openPMD::RecordComponent::SCALAR];
                    shadowgram.resetDataset(dataset);

                    // do not delete this object before dataPtr is not required anymore
                    auto data = helper->getShadowgramBuf();
                    auto sharedDataPtr = std::shared_ptr<float_64>{data->getPointer(), [](auto const*) {}};

                    shadowgram.storeChunk(sharedDataPtr, offset, extent);

                    series.iterations[currentStep].close();
                }


                void writeFile(std::vector<std::vector<float_64>> values, std::string name)
                {
                    std::ofstream outFile;
                    outFile.open(name.c_str(), std::ofstream::out | std::ostream::trunc);

                    if(!outFile)
                    {
                        std::cerr << "Can't open file [" << name << "] for output, disable plugin output. "
                                  << std::endl;
                    }
                    else
                    {
                        for(unsigned int i = 0; i < helper->getSizeX(); ++i) // over all x
                        {
                            for(unsigned int j = 0; j < helper->getSizeY(); ++j) // over all y
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

        } // namespace shadowgraphy
    } // namespace plugins
} // namespace picongpu