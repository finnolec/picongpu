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
#include "picongpu/plugins/shadowgraphy/ShadowgraphyHelper.hpp"

#include <pmacc/dataManagement/DataConnector.hpp>
#include <pmacc/math/Vector.hpp>
#include <pmacc/math/vector/Float.hpp>
#include <pmacc/math/vector/Int.hpp>
#include <pmacc/math/vector/Size_t.hpp>
#include <pmacc/mpi/MPIReduce.hpp>
#include <pmacc/mpi/reduceMethods/Reduce.hpp>

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

                MPI_Comm gatherComm = MPI_COMM_NULL;
                // gather rank zero is the master to dump files
                int gatherRank = -1;
                int numDevicesInPlane = 0;

                DataSpace<DIM2> localSliceSize;
                DataSpace<DIM2> localSliceOffset;
                DataSpace<DIM2> globalDomainSliceSize;
                int globalPlaneIdx = -1;
                int localPlaneIdx = -1;

                shadowgraphy::Helper* helper = nullptr;

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


                void buildGatherCommunicator(bool isActive)
                {
                    int countRanks;
                    int mpiRank;
                    MPI_CHECK(MPI_Comm_size(MPI_COMM_WORLD, &countRanks));
                    std::vector<int> allRank(countRanks);
                    std::vector<int> groupRanks(countRanks);
                    MPI_CHECK(MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank));

                    if(!isActive)
                        mpiRank = -1;

                    // avoid deadlock between not finished pmacc tasks and mpi blocking collectives
                    __getTransactionEvent().waitForFinished();
                    MPI_CHECK(MPI_Allgather(&mpiRank, 1, MPI_INT, allRank.data(), 1, MPI_INT, MPI_COMM_WORLD));

                    int numRanks = 0;
                    for(int i = 0; i < countRanks; ++i)
                    {
                        if(allRank[i] != -1)
                        {
                            std::cout << "allRank[i] i=" << i << " value=" << allRank[i] << std::endl;
                            groupRanks[numRanks] = allRank[i];
                            numRanks++;
                        }
                    }
                    numDevicesInPlane = numRanks;

                    MPI_Group group = MPI_GROUP_NULL;
                    MPI_Group newgroup = MPI_GROUP_NULL;
                    MPI_CHECK(MPI_Comm_group(MPI_COMM_WORLD, &group));
                    MPI_CHECK(MPI_Group_incl(group, numRanks, groupRanks.data(), &newgroup));

                    MPI_CHECK(MPI_Comm_create(MPI_COMM_WORLD, newgroup, &gatherComm));

                    if(mpiRank != -1)
                    {
                        MPI_CHECK(MPI_Comm_rank(gatherComm, &gatherRank));
                        std::cout << "gather rank=" << gatherRank << std::endl;
                    }
                    MPI_CHECK(MPI_Group_free(&group));
                    MPI_CHECK(MPI_Group_free(&newgroup));
                }

                void initGather()
                {
                    bool activatePlugin = true;

                    if(activatePlugin)
                    {
                        const SubGrid<simDim>& subGrid = Environment<simDim>::get().SubGrid();

                        auto globalDomain = subGrid.getGlobalDomain();
                        globalDomainSliceSize = DataSpace<DIM2>(globalDomain.size.x(), globalDomain.size.y());

                        auto globalPlaneExtent = globalDomain.size[plane];

                        globalPlaneIdx = globalPlaneExtent * slicePoint;

                        auto localDomain = subGrid.getLocalDomain();

                        auto isPlaneInLocalDomain = globalPlaneIdx >= localDomain.offset[plane]
                            && globalPlaneIdx < localDomain.offset[plane] + localDomain.size[plane];


                        if(isPlaneInLocalDomain)
                            localPlaneIdx = globalPlaneIdx - localDomain.offset[plane];

                        std::cout << "isPlaneInLocalDomain=" << isPlaneInLocalDomain
                                  << " localPlaneIdx=" << localPlaneIdx << " globalPlaneExtent=" << globalPlaneExtent
                                  << std::endl;


                        // collective call
                        buildGatherCommunicator(isPlaneInLocalDomain);
                        localSliceSize = DataSpace<DIM2>(localDomain.size.x(), localDomain.size.y());
                        localSliceOffset = DataSpace<DIM2>(localDomain.offset.x(), localDomain.offset.y());
                        std::cout << "localSliceOffset=" << localSliceOffset.toString() << std::endl;
                        std::cout << "localSliceSize=" << localSliceSize.toString() << std::endl;
                    }
                }

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

                            initGather();
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
                    if(gatherRank == -1)
                        return;
                    /* notification callback for simulation step currentStep
                     * called every notifyPeriod steps */
                    if(sliceIsOK)
                    {
                        bool isMaster = gatherRank == 0;

                        // First time the plugin is called:
                        if(isIntegrating == false)
                        {
                            startTime = currentStep;

                            if(isMaster)
                            {
                                std::cout << "master init " << currentStep << std::endl;
                                /// @todo shared pointer
                                helper = new Helper(
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
                            std::cout << "calculate " << currentStep << std::endl;

                            DataConnector& dc = Environment<>::get().DataConnector();
                            auto fieldE = dc.get<FieldE>(FieldE::getName(), false);

                            std::cout << "prepare E field" << std::endl;

                            storeSlice<shadowgraphy::Helper::FieldType::E>(
                                fieldE->getHostDataBox().shift(m_cellDescription->getGridLayout().getGuard()),
                                this->plane,
                                this->slicePoint,
                                localStep,
                                currentStep);

                            auto fieldB = dc.get<FieldB>(FieldB::getName(), false);

                            std::cout << "prepare B field" << std::endl;

                            storeSlice<shadowgraphy::Helper::FieldType::B>(
                                fieldB->getHostDataBox().shift(m_cellDescription->getGridLayout().getGuard()),
                                this->plane,
                                this->slicePoint,
                                localStep,
                                currentStep);


                            if(isMaster)
                            {
                                helper->calculate_dft(localStep);
                            }
                        }
                        else
                        {
                            if(isMaster)
                            {
                                std::cout << "dump " << currentStep << std::endl;
                                helper->propagateFields();
                                helper->calculate_shadowgram();

                                std::ostringstream filename;
                                filename << this->fileName << "_" << startTime << ":" << currentStep << ".dat";

                                writeFile(helper->getShadowgram(), filename.str());

                                writeToOpenPMDFile(currentStep);

                                delete(helper);
                            }
                            isIntegrating = false;
                        }
                    }
                }

                /* Stores the field slices from the host on the device. 2 field slices are required to adjust for
                 * the Yee-offset in the plugin.
                 * https://picongpu.readthedocs.io/en/latest/models/AOFDTD.html#maxwell-s-equations-on-the-mesh
                 * The current implementation is based on the (outdated) slice field printer printer and uses
                 * custl! It works, but it's not nice.
                 */
                template<typename shadowgraphy::Helper::FieldType T_fieldType, typename T_FieldBox>
                void storeSlice(const T_FieldBox& field, int nAxis, float slicePoint, int localStep, int currentStep)
                {
                    bool isMaster = gatherRank == 0;

                    auto pointField = std::make_shared<HostBufferIntern<float2_X, DIM2>>(localSliceSize);
                    auto pointFieldBox = pointField->getDataBox();

                    std::cout << "[" << gatherRank << "]"
                              << " start loading slice" << std::endl;
                    for(int y = 0; y < localSliceSize.y(); ++y)
                        for(int x = 0; x < localSliceSize.x(); ++x)
                        {
                            DataSpace<DIM2> idx(x, y);
                            DataSpace<DIM3> srcIdx(idx.x(), idx.y(), localPlaneIdx);
                            pointFieldBox(idx) = helper->cross<T_fieldType>(field.shift(srcIdx));
                        }
                    std::cout << "[" << gatherRank << "]"
                              << " end loding slice" << std::endl;

                    pmacc::GridController<simDim>& con = pmacc::Environment<simDim>::get().GridController();
                    auto numDevices = con.getGpuNodes();

                    std::cout << "[" << gatherRank << "]"
                              << " num devices in slice plane =" << numDevicesInPlane << std::endl;

                    // avoid deadlock between not finished pmacc tasks and mpi blocking collectives
                    __getTransactionEvent().waitForFinished();
                    // get number of elements per participating mpi rank
                    auto extentPerDevice = std::vector<DataSpace<DIM2>>(numDevicesInPlane);

                    std::cout << "[" << gatherRank << "]"
                              << " start gather extents " << localSliceSize.toString() << std::endl;
                    // gather extents
                    MPI_CHECK(MPI_Gather(
                        reinterpret_cast<int*>(&localSliceSize),
                        2,
                        MPI_INT,
                        reinterpret_cast<int*>(extentPerDevice.data()),
                        2,
                        MPI_INT,
                        0,
                        gatherComm));

                    if(isMaster)
                    {
                        for(int i = 0; i < numDevicesInPlane; ++i)
                        {
                            std::cout << "[" << gatherRank << "]"
                                      << " extent recive=" << extentPerDevice[i].toString() << std::endl;
                        }
                    }

                    std::cout << "[" << gatherRank << "]"
                              << " end gather extents" << std::endl;

                    auto offsetPerDevice = std::vector<DataSpace<DIM2>>(numDevicesInPlane);

                    std::cout << "[" << gatherRank << "]"
                              << " start gather offsets " << localSliceOffset.toString() << std::endl;

                    // gather offsets
                    MPI_CHECK(MPI_Gather(
                        reinterpret_cast<int*>(&localSliceOffset),
                        2,
                        MPI_INT,
                        reinterpret_cast<int*>(offsetPerDevice.data()),
                        2,
                        MPI_INT,
                        0,
                        gatherComm));

                    if(isMaster)
                    {
                        for(int i = 0; i < numDevicesInPlane; ++i)
                        {
                            std::cout << "[" << gatherRank << "]"
                                      << " offset recive=" << offsetPerDevice[i].toString() << std::endl;
                        }
                    }

                    std::cout << "[" << gatherRank << "]"
                              << " end gather offsets" << std::endl;

                    std::vector<int> displs(numDevicesInPlane);
                    std::vector<int> count(numDevicesInPlane);
                    // @todo replace by std::scan
                    int offset = 0;
                    int globalNumElements = 0u;

                    if(isMaster)
                    {
                        for(int i = 0; i < numDevicesInPlane; ++i)
                        {
                            std::cout << "[" << gatherRank << "]"
                                      << " offset=" << offset << std::endl;

                            displs[i] = offset * sizeof(float2_X);
                            count[i] = extentPerDevice[i].productOfComponents() * sizeof(float2_X);
                            offset += extentPerDevice[i].productOfComponents();
                            globalNumElements += extentPerDevice[i].productOfComponents();

                            std::cout << "[" << gatherRank << "]"
                                      << " extentPerDevice[" << i << "]=" << extentPerDevice[i] << std::endl;
                            std::cout << "[" << gatherRank << "]"
                                      << " displs[" << i << "]=" << displs[i] << std::endl;
                            std::cout << "[" << gatherRank << "]"
                                      << " count[" << i << "]=" << count[i] << std::endl;
                        }
                    }
                    std::cout << "[" << gatherRank << "]"
                              << " globalNumElements=" << globalNumElements << std::endl;

                    // gather all data from other ranks
                    auto allData = std::vector<float2_X>(globalNumElements);
                    int localNumElements = localSliceSize.productOfComponents();

                    MPI_CHECK(MPI_Gatherv(
                        reinterpret_cast<char*>(pointFieldBox.getPointer()),
                        localNumElements * sizeof(float2_X),
                        MPI_CHAR,
                        reinterpret_cast<char*>(allData.data()),
                        count.data(),
                        displs.data(),
                        MPI_CHAR,
                        0,
                        gatherComm));

                    std::cout << "[" << gatherRank << "]"
                              << " finish MPI_Gatherv" << std::endl;


                    if(isMaster)
                    {
                        auto globalField = std::make_shared<HostBufferIntern<float2_X, DIM2>>(globalDomainSliceSize);
                        auto globalFieldBox = globalField->getDataBox();

                        // aggregate data of all MPI ranks into a single 2D buffer
                        for(int dataSetNumber = 0; dataSetNumber < numDevicesInPlane; ++dataSetNumber)
                        {
                            for(int y = 0; y < extentPerDevice[dataSetNumber].y(); ++y)
                                for(int x = 0; x < extentPerDevice[dataSetNumber].x(); ++x)
                                {
                                    globalFieldBox(DataSpace<DIM2>(x, y) + offsetPerDevice[dataSetNumber]) = allData
                                        [displs[dataSetNumber] / sizeof(float2_X)
                                         + y * extentPerDevice[dataSetNumber].x() + x];
                                }
                        }

                        std::cout << "[" << gatherRank << "]"
                                  << " finish preparing global slice" << std::endl;
                        helper->storeField<T_fieldType>(localStep, currentStep, globalFieldBox);
                        std::cout << "[" << gatherRank << "]"
                                  << " finish store field" << std::endl;
                    }
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