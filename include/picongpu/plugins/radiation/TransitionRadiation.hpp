#pragma once

#include "picongpu/simulation_defines.hpp"

#include "picongpu/plugins/radiation/Radiation.kernel" // TODO Replace with TR Kernel
#include "picongpu/plugins/ISimulationPlugin.hpp"


#include <pmacc/mpi/reduceMethods/Reduce.hpp>
#include <pmacc/mpi/MPIReduce.hpp>
#include <pmacc/dimensions/DataSpaceOperations.hpp>
#include <pmacc/mappings/kernel/AreaMapping.hpp>

#include <boost/filesystem.hpp>

#include <string>
#include <iostream>
#include <cstdlib>

namespace picongpu
{
using namespace pmacc;

namespace po = boost::program_options;


///////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////  Transition Radiation Plugin Class  ///////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////
template<class ParticlesType>
class TransitionRadiation : public ISimulationPlugin
{
private:
  
    typedef MappingDesc::SuperCellSize SuperCellSize;

    typedef PIConGPUVerboseRadiation radLog;

    MappingDesc *cellDescription;
    std::string notifyPeriod;

    std::string speciesName;
    std::string pluginName;
    std::string pluginPrefix;
    std::string folderTotalRad;

    Amplitude* timeSumArray;

    bool isMaster;

    mpi::MPIReduce reduce;

public:
    // Constructor
    TransitionRadiation() :
    pluginName("TransitionRadiation: calculate transition radiation of species"),
    speciesName(ParticlesType::FrameType::getName()),
    pluginPrefix(speciesName + std::string("_transRad")),
    cellDescription(nullptr),
    timeSumArray(nullptr),
    isMaster(false)
    {
        log<radLog::SIMULATION_STATE> ("TR object created\n");

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
        log<radLog::SIMULATION_STATE> ("Notify called\n");
    }

    void pluginRegisterHelp(po::options_description& desc)
    {
        log<radLog::SIMULATION_STATE>("PluginRegisterHelp called\n");
        desc.add_options()
        ((pluginPrefix + ".period").c_str(), po::value<std::string> (&notifyPeriod), "enable plugin [for each n-th step]")
        ((pluginPrefix + ".folderTotalRad").c_str(), po::value<std::string > (&folderTotalRad)->default_value("totalRad"), "folder in which the integrated radiation from start of simulation is written");
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

    }

    void checkpoint(uint32_t timeStep, const std::string restartDirectory)
    {

    }

private:
    void pluginLoad()
    {
        log<radLog::SIMULATION_STATE>("Plugin Loaded");
        if(!notifyPeriod.empty())
        {
            Environment<>::get().PluginConnector().setNotificationPeriod(this, notifyPeriod);

            /*only rank 0 create a file*/
            isMaster = reduce.hasResult(mpi::reduceMethods::Reduce());
            pmacc::Filesystem<simDim>& fs = Environment<simDim>::get().Filesystem();

            if(isMaster)
            {
                fs.createDirectory("transitionRadiationHDF5");
                fs.setDirectoryPermissions("transitionRadiationHDF5");
            }


            if (isMaster)
            {
                //create folder for total output
                fs.createDirectory(folderTotalRad);
                fs.setDirectoryPermissions(folderTotalRad);
                for (unsigned int i = 0; i < elements_amplitude(); ++i)
                    timeSumArray[i] = Amplitude::zero();
            }
        }
    }

    void pluginUnload()
    {

    }

    void copyRadiationDeviceToHost()
    {

    }

    void saveRadPerGPU(const DataSpace<simDim> currentGPUpos)
    {

    }

    static unsigned int elements_amplitude()
    {
        return 0;
    }

    void collectRadiationOnMaster()
    {

    }

    void sumAmplitudesOverTime(Amplitude* targetArray, Amplitude* summandArray)
    {

    }

    void writeLastRadToText()
    {

    }

    void writeTotalRadToText()
    {

    }

    /** write total radiation data as HDF5 file */
    void writeAmplitudesToHDF5()
    {
        if (isMaster)
        {
            writeHDF5file(timeSumArray, std::string("transRadHDF5/") + speciesName + std::string("_radAmplitudes_"));
        }
    }

    void collectDataGPUToMaster()
    {

    }

    void writeAllFiles(const DataSpace<simDim> currentGPUpos)
    {

    }

    static const std::string dataLabels(int index)
    {
        return "0";
    }

    static const std::string dataLabelsDetectorDirection(int index)
    {
        return "0";
    }

    static const std::string dataLabelsDetectorFrequency(int index)
    {
        return "0";
    }

    static const std::string meshRecordLabels(int index)
    {
        return "0";
    }

    void writeHDF5file(Amplitude* values, std::string name)
    {

    }

    void readHDF5file(Amplitude* values, std::string name, const int timeStep)
    {

    }

    void writeFile(Amplitude* values, std::string name)
    {

    }

    void calculateRadiationParticles(uint32_t currentStep)
    {

    }

};

} // namespace picongpu