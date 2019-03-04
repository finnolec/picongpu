#pragma once

#include "picongpu/simulation_defines.hpp"

#include "picongpu/plugins/radiation/Radiation.kernel" // TODO Replace with TR Kernel
#include "picongpu/plugins/ISimulationPlugin.hpp"

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

public:
  // Constructor
  TransitionRadiation() :
  pluginName("TransitionRadiation: calculate transition radiation of species"),
  speciesName(ParticlesType::FrameType::getName()),
  pluginPrefix(speciesName + std::string("_transrad")),
  cellDescription(nullptr)
  {
    std::cout<<"TR object created";
    
    Environment<>::get().PluginConnector().registerPlugin(this);
    Environment<>::get().PluginConnector().setNotificationPeriod(this, "100");
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
    log<radLog::SIMULATION_STATE> ("Hello World, shoot\n");
    std::cout<<"Hello World\n";
  }

  void pluginRegisterHelp(po::options_description& desc)
  {
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