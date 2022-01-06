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
            }

            void setMappingDescription(MappingDesc *cellDescription) override
            {
            }

            private:
            std::string pluginName;
            std::string pluginPrefix;
            
            MappingDesc* cellDescription;
            std::string notifyPeriod;

            void pluginLoad() override
            {
                /* called when plugin is loaded, command line flags are available here
                * set notification period for our plugin at the PluginConnector */
                Environment<>::get().PluginConnector().setNotificationPeriod(this, this->notifyPeriod);
            }

            void pluginUnload() override
            {
                /* called when plugin is unloaded, cleanup here */
            }
        };
    }
}