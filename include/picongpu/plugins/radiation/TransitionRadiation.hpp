#pragma once

#include "picongpu/simulation_defines.hpp"

#include "picongpu/plugins/radiation/TransitionRadiation.kernel" // TODO Replace with TR Kernel
#include "picongpu/plugins/ISimulationPlugin.hpp"
#include "picongpu/plugins/common/stringHelpers.hpp"

#include <pmacc/mpi/reduceMethods/Reduce.hpp>
#include <pmacc/mpi/MPIReduce.hpp>
#include <pmacc/nvidia/functors/Add.hpp>
#include <pmacc/dimensions/DataSpaceOperations.hpp>
#include <pmacc/dataManagement/DataConnector.hpp>
#include <pmacc/mappings/kernel/AreaMapping.hpp>
#include <pmacc/traits/HasIdentifier.hpp>
#include <pmacc/traits/GetNumWorkers.hpp>

#include <boost/filesystem.hpp>

#include <string>
#include <iostream>
#include <fstream>
#include <cstdlib>

namespace picongpu
{
using namespace pmacc;

namespace po = boost::program_options;

///////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////~///  Transition Radiation Plugin Class  ///~///////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////
template<class ParticlesType>
class TransitionRadiation : public ISimulationPlugin
{
private:
  
    typedef MappingDesc::SuperCellSize SuperCellSize;

    typedef PIConGPUVerboseRadiation radLog;

    GridBuffer<picongpu::float_X, DIM1> *incTransRad;
    GridBuffer<picongpu::float_X, DIM1> *numParticles;

    radiation_frequencies::InitFreqFunctor freqInit;
    radiation_frequencies::FreqFunctor freqFkt;

    picongpu::float_X *tmp_itr;
    picongpu::float_X *tmp_num;
    picongpu::float_X* theTransRad; 
    MappingDesc *cellDescription;
    std::string notifyPeriod;
    uint32_t timeStep;
    uint32_t potatoSalad;

    std::string speciesName;
    std::string pluginName;
    std::string pluginPrefix;
    std::string filename_prefix;
    std::string folderTransRad;
    std::string pathOmegaList;

    vector_X* detectorPositions;
    float_X* detectorFrequencies;

    bool isMaster;
    uint32_t currentStep;


    mpi::MPIReduce reduce;

    bool debug;

public:
    // Constructor
    TransitionRadiation() :
    pluginName("TransitionRadiation: calculate transition radiation of species"),
    speciesName(ParticlesType::FrameType::getName()),
    pluginPrefix(speciesName + std::string("_transRad")),
    folderTransRad("transRad"),
    filename_prefix(pluginPrefix),
    incTransRad(nullptr),
    numParticles(nullptr),
    cellDescription(nullptr),
    tmp_itr(nullptr),
    tmp_num(nullptr),
    theTransRad(nullptr),
    detectorPositions(nullptr),
    detectorFrequencies(nullptr),
    isMaster(false),
    currentStep(0),
    potatoSalad(100),
    debug(true)
    {
        if (debug)
        {
            log<radLog::SIMULATION_STATE> ("TR object created");
        }

        Environment<>::get().PluginConnector().registerPlugin(this);
        //Environment<>::get().PluginConnector().setNotificationPeriod(this, "100");
    }

    virtual ~TransitionRadiation()
    {
    }

    /**
     * This function calculates the Transition Radiation by calling the
     * according function of the Kernel File.
     * @param currentStep
     */
    void notify(uint32_t currentStep)
    {
        //log<radLog::SIMULATION_STATE> ("Transition Radition (%1%): calculate time step %2% ") % speciesName % currentStep;
        if (currentStep == potatoSalad)
        {
            //log<radLog::SIMULATION_STATE> ("Transition Radition (%1%): calculate time step %2% ") % speciesName % currentStep;
            if (debug)
            {
                log<radLog::SIMULATION_STATE> ("Notify called in if statement.");
            }
            calculateRadiationParticles(currentStep);
            if (debug)
            {
                log<radLog::SIMULATION_STATE> ("Notify finished in if statement.");
            }
        }
    }

    void pluginRegisterHelp(po::options_description& desc)
    {
        desc.add_options()
            ((pluginPrefix + ".period").c_str(), po::value<std::string> (&notifyPeriod), "enable plugin [for each n-th step]")
            ((pluginPrefix + ".omegaList").c_str(), po::value<std::string > (&pathOmegaList)->default_value("_noPath_"), "path to file containing all frequencies to calculate");
    // log<radLog::SIMULATION_STATE>("Notify period is: (%1%)" % &notifyPeriod);
    }

    std::string pluginGetName() const
    {
        return pluginName;
    }

    void setMappingDescription(MappingDesc *cellDescription)
    {
        this->cellDescription = cellDescription;
    }

    void restart(uint32_t timeStep, const std::string restartDirectory)
    {
        if (debug)
        {
            log<radLog::SIMULATION_STATE> ("restart called");
        }
    }

    void checkpoint(uint32_t timeStep, const std::string restartDirectory)
    {
        if (debug)
        {
            log<radLog::SIMULATION_STATE> ("checkpoint called");
        }
    }

private:
    void pluginLoad()
    {
        if(!notifyPeriod.empty())
        {
            log<radLog::SIMULATION_STATE>("Plugin Loaded");

            /*only rank 0 create a file*/
            isMaster = reduce.hasResult(mpi::reduceMethods::Reduce());
            pmacc::Filesystem<simDim>& fs = Environment<simDim>::get().Filesystem();

            tmp_itr = new picongpu::float_X[elements_amplitude()];
            tmp_num = new picongpu::float_X[elements_amplitude()];
            
            Environment<>::get().PluginConnector().setNotificationPeriod(this, notifyPeriod);

            incTransRad = new GridBuffer<picongpu::float_X, DIM1 > (DataSpace<DIM1 > (elements_amplitude()));
            numParticles = new GridBuffer<picongpu::float_X, DIM1 > (DataSpace<DIM1 > (elements_amplitude()));

            freqInit.Init(pathOmegaList);
            freqFkt = freqInit.getFunctor();

            if (isMaster)
            {
                theTransRad = new picongpu::float_X[elements_amplitude()];
                /* save detector position / observation direction */
                detectorPositions = new vector_X[parameters::N_observer];
                for(uint32_t detectorIndex=0; detectorIndex < parameters::N_observer; ++detectorIndex)
                {
                    detectorPositions[detectorIndex] = radiation_observer::observation_direction(detectorIndex);
                }

                /* save detector frequencies */
                detectorFrequencies = new float_X[radiation_frequencies::N_omega];
                for(uint32_t detectorIndex=0; detectorIndex < radiation_frequencies::N_omega; ++detectorIndex)
                {
                    detectorFrequencies[detectorIndex] = freqFkt(detectorIndex);
                }

                for (unsigned int i=0; i< elements_amplitude(); ++i)
                {
                    theTransRad[i] = 0;
                }

                fs.createDirectory(folderTransRad);
                fs.setDirectoryPermissions(folderTransRad);
            }
        }
    }

    void pluginUnload()
    {
        if (debug)
        {
            log<radLog::SIMULATION_STATE> ("pluginunload called");
        }
        collectDataGPUToMaster();
        writeTransRadToText();
        if(isMaster)
        {
            __deleteArray(theTransRad);
        }
        CUDA_CHECK(cudaGetLastError());
        __delete( incTransRad );
        __delete( numParticles );
        __deleteArray(tmp_itr);
        __deleteArray(tmp_num);
    }

    void copyRadiationDeviceToHost()
    {
        if(debug)
        {
            log<radLog::SIMULATION_STATE>("copyradiationdevicetohost called");
        }
        incTransRad->deviceToHost();
        numParticles->deviceToHost();
        __getTransactionEvent().waitForFinished();
    }

    static unsigned int elements_amplitude()
    {
        return radiation_frequencies::N_omega * parameters::N_observer; // storage for amplitude results on GPU
    }

    /** combine radiation data from each CPU and store result on master
     *  copyRadiationDeviceToHost() should be called before */
    void collectRadiationOnMaster()
    {
        if(debug)
        {
            log<radLog::SIMULATION_STATE>("collectRadiationOnMaster");
        }
        reduce(
            nvidia::functors::Add(),
            tmp_itr,
            incTransRad->getHostBuffer().getBasePointer(),
            elements_amplitude(),
            mpi::reduceMethods::Reduce()
        );
        reduce(
            nvidia::functors::Add(),
            tmp_num,
            numParticles->getHostBuffer().getBasePointer(),
            elements_amplitude(),
            mpi::reduceMethods::Reduce()
        );
    }

    void writeTransRadToText()
    {
        if(debug)
        {
            log<radLog::SIMULATION_STATE>("writeTransRadToText");
        }
        // only the master rank writes data
        if (isMaster)
        {
            // get time step as string
            std::stringstream o_step;
            o_step << currentStep;

            // write totalRad data to txt
            writeFile(theTransRad, folderTransRad + "/" + filename_prefix + "_" + o_step.str() + ".dat");
        }
    }


    /** perform all operations to get data from GPU to master */
    void collectDataGPUToMaster()
    {
        // collect data GPU -> CPU -> Master
        copyRadiationDeviceToHost();
        collectRadiationOnMaster();
        //sumTransitionRadiation(theTransRad, incTransRad->getHostBuffer().getBasePointer(), numParticles->getHostBuffer().getBasePointer());
        sumTransitionRadiation(theTransRad, tmp_itr, tmp_num);
    }

    // calculate transition radiation integrals with the energy values from the kernel
    void sumTransitionRadiation(
        picongpu::float_X* targetArray, 
        picongpu::float_X* itrArray,
        picongpu::float_X* numArray
    )
    {
        if (isMaster)
        {
            for(unsigned int i = 0; i < elements_amplitude(); ++i)
            {
                targetArray[i] = itrArray[i]; // * numArray[i];
            }
        }
    }


    static const std::string dataLabels(int index)
    {
        log<radLog::SIMULATION_STATE>("dataLabels");
    
        return "0";
    }

    static const std::string dataLabelsDetectorDirection(int index)
    {
    
        log<radLog::SIMULATION_STATE>("dataLabelsDetectorDirection");
        
        return "0";
    }

    static const std::string dataLabelsDetectorFrequency(int index)
    {
        log<radLog::SIMULATION_STATE>("dataLabelsDetectorFrequency");
        
        const std::string path("DetectorFrequency/");

        /* return record name if handed -1 */
        if(index == -1)
            return path;

        const std::string dataLabelsList[] = {"omega"};

        return path + dataLabelsList[index];
    }

    void writeFile(float_X* values, std::string name)
    {

        if(debug)
        {
            log<radLog::SIMULATION_STATE>("writeFile");
        }
        std::ofstream outFile;
        outFile.open(name.c_str(), std::ofstream::out | std::ostream::trunc);
        if (!outFile)
        {
            std::cerr << "Can't open file [" << name << "] for output, disable plugin output. " << std::endl;
            isMaster = false; // no Master anymore -> no process is able to write
        }
        else
        {
            for (unsigned int index_direction = 0; index_direction < parameters::N_observer; ++index_direction) // over all directions
            {
                for (unsigned index_omega = 0; index_omega < radiation_frequencies::N_omega; ++index_omega) // over all frequencies
                {
                    // Take Amplitude for one direction and frequency,
                    // calculate the square of the absolute value
                    // and write to file.
                    constexpr picongpu::float_X transRadUnit = 
                        UNIT_CHARGE * UNIT_CHARGE * (1.0 / (4 * PI * SI::EPS0_SI * PI * PI * SI::SPEED_OF_LIGHT_SI));
                    outFile <<
                        values[index_direction * radiation_frequencies::N_omega + index_omega] 
                        * transRadUnit << "\t";

                }
                outFile << std::endl;
            }
            outFile.flush();
            outFile << std::endl; //now all data are written to file

            if (outFile.fail())
                std::cerr << "Error on flushing file [" << name << "]. " << std::endl;

            outFile.close();
        }
    }

    void calculateRadiationParticles(uint32_t currentStep)
    {

        if(debug)
        {
            log<radLog::SIMULATION_STATE>("calculateRadiationParticles called");
        }

        DataConnector &dc = Environment<>::get().DataConnector();
        auto particles = dc.get< ParticlesType >( ParticlesType::FrameType::getName(), true);

        const int N_observer = parameters::N_observer;
        const auto gridDim_rad = N_observer;

        /* number of threads per block = number of cells in a super cell
        *          = number of particles in a Frame
        *          (THIS IS PIConGPU SPECIFIC)
        * A Frame is the entity that stores particles.
        * A super cell can have many Frames.
        * Particles in a Frame can be accessed in parallel.
        */

        // Some funny things that make it possible for the kernel to calculate
        // the absolute position of the particles
        DataSpace<simDim> localSize(cellDescription->getGridLayout().getDataSpaceWithoutGuarding());
        const uint32_t numSlides = MovingWindow::getInstance().getSlideCounter(currentStep);
        const SubGrid<simDim>& subGrid = Environment<simDim>::get().SubGrid();
        DataSpace<simDim> globalOffset(subGrid.getLocalDomain().offset);
        globalOffset.y() += (localSize.y() * numSlides);

        constexpr uint32_t numWorkers = pmacc::traits::GetNumWorkers<
            pmacc::math::CT::volume< SuperCellSize >::type::value
        >::value;

        // PIC-like kernel call of the radiation kernel
        PMACC_KERNEL( KernelTransRadParticles<
            numWorkers
        >{} )(
            gridDim_rad,
            numWorkers
        )(
            /*Pointer to particles memory on the device*/
            particles->getDeviceParticlesBox(),

            /*Pointer to memory of radiated amplitude on the device*/
            incTransRad->getDeviceBuffer().getDataBox(),
            numParticles->getDeviceBuffer().getDataBox(),
            globalOffset,
            *cellDescription,
            freqFkt,
            subGrid.getGlobalDomain().size
        );

        dc.releaseData( ParticlesType::FrameType::getName() );

        if(debug)
        {
            log<radLog::SIMULATION_STATE>("calculateRadiationParticles finished");
        }

    }
};

} // namespace picongpu