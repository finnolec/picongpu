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

    GridBuffer<picongpu::float_64, DIM1> *radiation;
    GridBuffer<picongpu::float_64, DIM1> *incTransRad;
    GridBuffer<picongpu::float_64, DIM1> *cohTransRadPerp;
    GridBuffer<picongpu::float_64, DIM1> *cohTransRadPara;

    radiation_frequencies::InitFreqFunctor freqInit;
    radiation_frequencies::FreqFunctor freqFkt;

    MappingDesc *cellDescription;
    std::string notifyPeriod;

    std::string speciesName;
    std::string pluginName;
    std::string pluginPrefix;
    std::string filename_prefix;
    std::string folderTransRad;
    std::string pathOmegaList;

    vector_64* detectorPositions;
    float_64* detectorFrequencies;

    bool isMaster;

    bool debug;

public:
    // Constructor
    TransitionRadiation() :
    pluginName("TransitionRadiation: calculate transition radiation of species"),
    speciesName(ParticlesType::FrameType::getName()),
    pluginPrefix(speciesName + std::string("_transRad")),
    filename_prefix(pluginPrefix),
    radiation(nullptr),
    cellDescription(nullptr),
    totalRad(false),
    timeSumArray(nullptr),
    tmp_result(nullptr),
    detectorPositions(nullptr),
    detectorFrequencies(nullptr),
    isMaster(false),
    currentStep(0),
    lastStep(0),
    compressionOn(false),
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
     * @param: currentStep
     */
    void notify(uint32_t currentStep)
    {
        //log<radLog::SIMULATION_STATE> ("Transition Radition (%1%): calculate time step %2% ") % speciesName % currentStep;
        if (debug)
        {
            log<radLog::SIMULATION_STATE> ("Notify called");
        }
    }

    void pluginRegisterHelp(po::options_description& desc)
    {
        if (debug)
        {
            log<radLog::SIMULATION_STATE>("PluginRegisterHelp called");
        }
        desc.add_options()
        ((pluginPrefix + ".period").c_str(), po::value<std::string> (&notifyPeriod), "enable plugin [for each n-th step]")
        ((pluginPrefix + ".totalRadiation").c_str(), po::bool_switch(&totalRad), "enable calculation of integrated radiation from start of simulation")
        ((pluginPrefix + ".folderTotalRad").c_str(), po::value<std::string > (&folderTotalRad)->default_value("totalRad"), "folder in which the integrated radiation from start of simulation is written")
        ((pluginPrefix + ".omegaList").c_str(), po::value<std::string > (&pathOmegaList)->default_value("_noPath_"), "path to file containing all frequencies to calculate")
        ((pluginPrefix + ".compression").c_str(), po::bool_switch(&compressionOn), "enable compression of hdf5 output");

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
        log<radLog::SIMULATION_STATE>("Plugin Loaded");
        if(!notifyPeriod.empty())
        {
            Environment<>::get().PluginConnector().setNotificationPeriod(this, notifyPeriod);

            radiation = new GridBuffer<Amplitude, DIM1 > (DataSpace<DIM1 > (elements_amplitude())); //create one int on GPU and host

            freqInit.Init(pathOmegaList);
            freqFkt = freqInit.getFunctor();

            /*only rank 0 create a file*/
            isMaster = reduce.hasResult(mpi::reduceMethods::Reduce());
            pmacc::Filesystem<simDim>& fs = Environment<simDim>::get().Filesystem();

            if (isMaster)
            {
                timeSumArray = new Amplitude[elements_amplitude()];

                /* save detector position / observation direction */
                detectorPositions = new vector_64[parameters::N_observer];
                for(uint32_t detectorIndex=0; detectorIndex < parameters::N_observer; ++detectorIndex)
                {
                    detectorPositions[detectorIndex] = radiation_observer::observation_direction(detectorIndex);
                }

                /* save detector frequencies */
                detectorFrequencies = new float_64[radiation_frequencies::N_omega];
                for(uint32_t detectorIndex=0; detectorIndex < radiation_frequencies::N_omega; ++detectorIndex)
                {
                    detectorFrequencies[detectorIndex] = freqFkt(detectorIndex);
                }

                if(totalRad)
                {
                    fs.createDirectory("transitionRadiationHDF5");
                    fs.setDirectoryPermissions("transitionRadiationHDF5");

                    //create folder for total output
                    fs.createDirectory(folderTotalRad);
                    fs.setDirectoryPermissions(folderTotalRad);
                    for (unsigned int i = 0; i < elements_amplitude(); ++i)
                        timeSumArray[i] = Amplitude::zero();
                }
            }
        }
    }

    void pluginUnload()
    {
        if (debug)
        {
            log<radLog::SIMULATION_STATE> ("pluginunload called");
        }
    }

    void copyRadiationDeviceToHost()
    {
        if(debug)
        {
            log<radLog::SIMULATION_STATE>("copyradiationdevicetohost called");
        }
        radiation->deviceToHost();
        __getTransactionEvent().waitForFinished();
    }

    void saveRadPerGPU(const DataSpace<simDim> currentGPUpos)
    {
        if(debug)
        {
            log<radLog::SIMULATION_STATE>("saveRadPerGPU called");
        }

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
        reduce(nvidia::functors::Add(),
                tmp_result,
                radiation->getHostBuffer().getBasePointer(),
                elements_amplitude(),
                mpi::reduceMethods::Reduce()
                );
    }


    /** add collected radiation data to previously stored data
     *  should be called after collectRadiationOnMaster() */
    void sumAmplitudesOverTime(Amplitude* targetArray, Amplitude* summandArray)
    {
        if(debug)
        {
            log<radLog::SIMULATION_STATE>("sumAmplitudesOverTime");
        }
        if (isMaster)
        {
            // add last amplitudes to previous amplitudes
            for (unsigned int i = 0; i < elements_amplitude(); ++i)
            targetArray[i] += summandArray[i];
        }
    }

    void writeLastRadToText()
    {
        if(debug)
        {
            log<radLog::SIMULATION_STATE>("writeLastRadToText");
        }
    }

    void writeTotalRadToText()
    {
        if(debug)
        {
            log<radLog::SIMULATION_STATE>("writeTotalRadToText");
        }
        // only the master rank writes data
        if (isMaster)
        {
            // write file only if totalRad flag was selected
            if (totalRad)
            {
                // get time step as string
                std::stringstream o_step;
                o_step << currentStep;

                // write totalRad data to txt
                writeFile(timeSumArray, folderTotalRad + "/" + filename_prefix + "_" + o_step.str() + ".dat");
            }
        }
    }

    /** write total radiation data as HDF5 file */
    void writeAmplitudesToHDF5()
    {
        if(debug)
        {
            log<radLog::SIMULATION_STATE>("writeAmplitudesToHDF5");
        }
        if (isMaster)
        {
            writeHDF5file(timeSumArray, std::string("transRadHDF5/") + speciesName + std::string("_radAmplitudes_"));
        }
    }

    /** perform all operations to get data from GPU to master */
    void collectDataGPUToMaster()
    {

        if(debug)
        {
            log<radLog::SIMULATION_STATE>("collectDataGPUToMaster");
        }
        // collect data GPU -> CPU -> Master
        copyRadiationDeviceToHost();
        collectRadiationOnMaster();
        sumAmplitudesOverTime(timeSumArray, tmp_result);
    }

    /** write all possible/selected output */
    void writeAllFiles(const DataSpace<simDim> currentGPUpos)
    {
        if(debug)
        {
            log<radLog::SIMULATION_STATE>("writeAllFiles");
        }
        // write data to files
        saveRadPerGPU(currentGPUpos);
        writeLastRadToText();
        writeTotalRadToText();
        writeAmplitudesToHDF5();
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

    static const std::string meshRecordLabels(int index)
    {
        log<radLog::SIMULATION_STATE>("meshRecordLabels");
        
        if(index == idLabels::Amplitude)
            return dataLabels(-1);
        else if (index == idLabels::Detector)
            return dataLabelsDetectorDirection(-1);
        else if (index == idLabels::Frequency)
            return dataLabelsDetectorFrequency(-1);
        else
            return std::string("this-record-does-not-exist");
    }

    void writeHDF5file(Amplitude* values, std::string name)
    {  
    
        if(debug)
        {
            log<radLog::SIMULATION_STATE>("writeHDF5file");
        }

    }

    void readHDF5file(Amplitude* values, std::string name, const int timeStep)
    {
        if(debug)
        {
            log<radLog::SIMULATION_STATE>("readHDF5file");
        }

    }

    void writeFile(Amplitude* values, std::string name)
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
                    outFile <<
                        values[index_omega + index_direction * radiation_frequencies::N_omega].calc_radiation() * UNIT_ENERGY * UNIT_TIME << "\t";

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
            log<radLog::SIMULATION_STATE>("calculateRadiationParticles");
        }

        this->currentStep = currentStep;
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
        PMACC_KERNEL( KernelRadiationParticles<
            numWorkers,
        >{} )(
            gridDim_rad,
            numWorkers
        )(
            /*Pointer to particles memory on the device*/
            particles->getDeviceParticlesBox(),

            /*Pointer to memory of radiated amplitude on the device*/
            incTransRad->getDeviceBuffer().getDataBox(),
            cohTransRadPara->getDeviceBuffer().getDataBox(),
            cohTransRadPerp->getDeviceBuffer().getDataBox(),
            numParticles->getDeviceBuffer().getDataBox(),
            globalOffset,
            *cellDescription,
            freqFkt,
            subGrid.getGlobalDomain().size
        );

        dc.releaseData( ParticlesType::FrameType::getName() );
        

    }

};

} // namespace picongpu