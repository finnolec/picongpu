#pragma once

#include "picongpu/simulation_defines.hpp"

#include "picongpu/plugins/radiation/TransitionRadiation.kernel"
#include "picongpu/plugins/ISimulationPlugin.hpp"
#include "picongpu/plugins/radiation/ExecuteParticleFilter.hpp"
#include "picongpu/plugins/common/stringHelpers.hpp"

#include <pmacc/mpi/reduceMethods/Reduce.hpp>
#include <pmacc/mpi/MPIReduce.hpp>
#include <pmacc/nvidia/functors/Add.hpp>
#include <pmacc/dimensions/DataSpaceOperations.hpp>
#include <pmacc/dataManagement/DataConnector.hpp>
#include <pmacc/mappings/kernel/AreaMapping.hpp>
#include <pmacc/traits/HasIdentifier.hpp>
#include <pmacc/traits/GetNumWorkers.hpp>

#include <pmacc/math/Complex.hpp>

#include <boost/filesystem.hpp>

#include <string>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <ctime>


namespace picongpu
{
    using namespace pmacc;

    namespace po = boost::program_options;
    using complex_X = pmacc::math::Complex< float_X >;

    ///////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////~///  Transition Radiation Plugin Class  ///~///////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////
    template<
        typename T_ParticlesType
    >
    class TransitionRadiation : public ISimulationPlugin
    {
    private:
    
        using SuperCellSize = MappingDesc::SuperCellSize;

        using radLog = PIConGPUVerboseRadiation;

        GridBuffer< float_X, DIM1 > * incTransRad;
        GridBuffer< complex_X, DIM1 > * cohTransRadPara;
        GridBuffer< complex_X, DIM1 > * cohTransRadPerp;
        GridBuffer< float_X, DIM1 > * numParticles;

        radiation_frequencies::InitFreqFunctor freqInit;
        radiation_frequencies::FreqFunctor freqFkt;

        float_X * tmpITR;
        complex_X * tmpCTRpara;
        complex_X * tmpCTRperp;
        float_X * tmpNum;
        float_X * theTransRad; 
        MappingDesc * cellDescription;
        std::string notifyPeriod;
        uint32_t timeStep;

        std::string speciesName;
        std::string pluginName;
        std::string pluginPrefix;
        std::string filenamePrefix;
        std::string folderTransRad;
        std::string pathOmegaList;

        float3_X * detectorPositions;
        float_X * detectorFrequencies;

        bool isMaster;
        uint32_t currentStep;

        mpi::MPIReduce reduce;

    public:
        // Constructor
        TransitionRadiation( ) :
            pluginName( "TransitionRadiation: calculate transition radiation of species" ),
            speciesName( T_ParticlesType::FrameType::getName( ) ),
            pluginPrefix( speciesName + std::string( "_transRad" ) ),
            folderTransRad( "transRad" ),
            filenamePrefix( pluginPrefix ),
            incTransRad( nullptr ),
            cohTransRadPara( nullptr ),
            cohTransRadPerp( nullptr ),
            numParticles( nullptr ),
            cellDescription( nullptr ),
            tmpITR( nullptr ),
            tmpCTRpara( nullptr ),
            tmpCTRperp( nullptr ),
            tmpNum( nullptr ),
            theTransRad( nullptr ),
            detectorPositions( nullptr ),
            detectorFrequencies( nullptr ),
            isMaster( false ),
            currentStep( 0 )
        {
            Environment< >::get( ).PluginConnector( ).registerPlugin( this );
        }

        virtual 
        ~TransitionRadiation( )
        {
        }

        /**
         * This function calculates the Transition Radiation by calling the
         * according function of the Kernel File.
         * @param currentStep
         */
        void 
        notify(
            uint32_t currentStep
        )
        {
            log< radLog::SIMULATION_STATE >( "Transition Radition (%1%): calculate time step %2% " ) % speciesName % currentStep;
            
            resetBuffers( );
            this->currentStep = currentStep;
            
            // const std::clock_t beginTime = std::clock( );

            calculateRadiationParticles( currentStep );
            
            // std::cout << float( std::clock( ) - beginTime ) / std::CLOCKS_PER_SEC;

            log< radLog::SIMULATION_STATE >( "Transition Radition (%1%): finished time step %2% " ) % speciesName % currentStep;
            
            collectDataGPUToMaster( );
            writeTransRadToText( );

            // std::cout << float( std::clock( ) - beginTime ) / std::CLOCKS_PER_SEC;

            
            log< radLog::SIMULATION_STATE >( "Transition Radition (%1%): printed to table %2% " ) % speciesName % currentStep;
        }

        void 
        pluginRegisterHelp(
            po::options_description& desc
        )
        {
            desc.add_options( )(
                ( pluginPrefix + ".period" ).c_str( ), 
                po::value<std::string> ( &notifyPeriod ), 
                "enable plugin [for each n-th step]"
            );
        }

        std::string 
        pluginGetName( ) const
        {
            return pluginName;
        }

        void 
        setMappingDescription(
            MappingDesc *cellDescription
        )
        {
            this->cellDescription = cellDescription;
        }

        void 
        restart(
            uint32_t timeStep, 
            const std::string restartDirectory
        )
        {
        }

        void 
        checkpoint(
            uint32_t timeStep, 
            const std::string restartDirectory
        )
        {
        }

    private:
        void 
        resetBuffers ( )
        {
            /* Resets all Databuffers and arrays for repeated calculation of the 
            * transition radiation
            */
            incTransRad->getDeviceBuffer( ).reset( false );
            cohTransRadPara->getDeviceBuffer( ).reset( false );
            cohTransRadPerp->getDeviceBuffer( ).reset( false );
            numParticles->getDeviceBuffer( ).reset( false );

            for( unsigned int i=0; i< elements_amplitude( ); ++i )
            {
                tmpITR[ i ] = 0;
                tmpCTRpara[ i ] = 0;
                tmpCTRperp[ i ] = 0;
                tmpNum[ i ] = 0;
                if( isMaster )
                {
                    theTransRad[ i ] = 0;  
                }
            }
        }

        void 
        pluginLoad( )
        {
            
            tmpITR = new float_X[ elements_amplitude( ) ];
            tmpCTRpara = new complex_X[ elements_amplitude( ) ];
            tmpCTRperp = new complex_X[ elements_amplitude( ) ];
            tmpNum = new float_X[ elements_amplitude( ) ];
            if( !notifyPeriod.empty( ) )
            {
                /*only rank 0 create a file*/
                isMaster = reduce.hasResult( mpi::reduceMethods::Reduce( ) );
                pmacc::Filesystem<simDim>& fs = Environment<simDim>::get( ).Filesystem( );

                Environment<>::get( ).PluginConnector( ).setNotificationPeriod( this, notifyPeriod );

                incTransRad = new GridBuffer< float_X, DIM1 >( 
                    DataSpace< DIM1 > ( elements_amplitude( ) ) );
                cohTransRadPara = new GridBuffer< complex_X, DIM1 >( 
                    DataSpace< DIM1 > ( elements_amplitude( ) ) );
                cohTransRadPerp = new GridBuffer< complex_X, DIM1 >( 
                    DataSpace< DIM1 > ( elements_amplitude( ) ) );
                numParticles = new GridBuffer< float_X, DIM1 >( 
                    DataSpace< DIM1 > ( elements_amplitude( ) ) );

                freqInit.Init( pathOmegaList );
                freqFkt = freqInit.getFunctor( );

                if ( isMaster )
                {
                    theTransRad = new float_X[ elements_amplitude( ) ];
                    /* save detector position / observation direction */
                    detectorPositions = new float3_X[ parameters::N_observer ];
                    for(
                        uint32_t detectorIndex=0; 
                        detectorIndex < parameters::N_observer; 
                        ++detectorIndex
                    )
                    {
                        detectorPositions[ detectorIndex ] = radiation_observer::observation_direction_picongpustandard( detectorIndex );
                    }

                    /* save detector frequencies */
                    detectorFrequencies = new float_X[ radiation_frequencies::N_omega ];
                    for(
                        uint32_t detectorIndex=0; 
                        detectorIndex < radiation_frequencies::N_omega; 
                        ++detectorIndex
                    )
                    {
                        detectorFrequencies[ detectorIndex ] = freqFkt( detectorIndex );
                    }

                    for ( unsigned int i=0; i< elements_amplitude( ); ++i )
                    {
                        theTransRad[ i ] = 0;
                    }

                    fs.createDirectory( folderTransRad );
                    fs.setDirectoryPermissions( folderTransRad );
                }
            }
        }

        void 
        pluginUnload( )
        {
            if( isMaster )
            {
                __deleteArray( theTransRad );
            }
            CUDA_CHECK( cudaGetLastError( ) );
            __delete( incTransRad );
            __delete( cohTransRadPara );
            __delete( cohTransRadPerp );
            __delete( numParticles );
            __deleteArray( tmpITR );
            __deleteArray( tmpCTRpara );
            __deleteArray( tmpCTRperp );
            __deleteArray( tmpNum );
        }

        void 
        copyRadiationDeviceToHost( )
        {
            incTransRad->deviceToHost( );
            __getTransactionEvent( ).waitForFinished( );
            cohTransRadPara->deviceToHost( );
            __getTransactionEvent( ).waitForFinished( );
            cohTransRadPerp->deviceToHost( );
            __getTransactionEvent( ).waitForFinished( );
            numParticles->deviceToHost( );
            __getTransactionEvent( ).waitForFinished( );
        }

        static 
        unsigned int 
        elements_amplitude( )
        {
            return radiation_frequencies::N_omega * parameters::N_observer; // storage for amplitude results on GPU
        }

        /** combine radiation data from each CPU and store result on master
         *  copyRadiationDeviceToHost( ) should be called before */
        void 
        collectRadiationOnMaster( )
        {
            reduce(
                nvidia::functors::Add( ),
                tmpITR,
                incTransRad->getHostBuffer( ).getBasePointer( ),
                elements_amplitude( ),
                mpi::reduceMethods::Reduce( )
            );
            reduce(
                nvidia::functors::Add( ),
                tmpCTRpara,
                cohTransRadPara->getHostBuffer( ).getBasePointer( ),
                elements_amplitude( ),
                mpi::reduceMethods::Reduce( )
            );
            reduce(
                nvidia::functors::Add( ),
                tmpCTRperp,
                cohTransRadPerp->getHostBuffer( ).getBasePointer( ),
                elements_amplitude( ),
                mpi::reduceMethods::Reduce( )
            );
            reduce(
                nvidia::functors::Add( ),
                tmpNum,
                numParticles->getHostBuffer( ).getBasePointer( ),
                elements_amplitude( ),
                mpi::reduceMethods::Reduce( )
            );
        }

        void 
        writeTransRadToText( )
        {
            // only the master rank writes data
            if (isMaster)
            {
                // get time step as string
                std::stringstream o_step;
                o_step << currentStep;

                // write totalRad data to txt
                writeFile(theTransRad, folderTransRad + "/" + filenamePrefix + "_" + o_step.str( ) + ".dat");
            }
        }


        /** perform all operations to get data from GPU to master */
        void 
        collectDataGPUToMaster( )
        {
            // collect data GPU -> CPU -> Master
            copyRadiationDeviceToHost( );
            collectRadiationOnMaster( );
            sumTransitionRadiation( theTransRad, tmpITR, tmpCTRpara, tmpCTRperp, tmpNum );
        }

        // calculate transition radiation integrals with the energy values from the kernel
        void 
        sumTransitionRadiation(
            float_X * targetArray, 
            float_X * itrArray,
            complex_X * ctrParaArray,
            complex_X * ctrPerpArray,
            float_X * numArray
        )
        {
            if (isMaster)
            {
                /************************************************************
                 ******** Here happens the true physical calculation ********
                ************************************************************/
                for( unsigned int i = 0; i < elements_amplitude( ); ++i )
                {
                    // std::cout << numArray[ i ] << " numArray[ i ]\n";
                    // std::cout << itrArray[ i ] << " itrArray[ i ]\n";
                    // std::cout << ctrParaArray[ i ].get_real( ) << " ctrParaArray[ i ]\n";
                    // std::cout << ctrPerpArray[ i ].get_imag( ) << " ctrPerpArray[ i ]\n";

                    const float_X ctrPara = math::abs2( ctrParaArray[ i ] );
                    const float_X ctrPerp = math::abs2( ctrPerpArray[ i ] );

                    targetArray[ i ] = ( 
                        itrArray[ i ] + ( numArray[ i ] - 1.0 ) * ( ctrPara + ctrPerp ) / numArray[i]
                    );
                    // targetArray[ i ] = ( 
                    //     itrArray[i]
                    // );
                }
            }
        }

        void 
        writeFile(
            float_X * values, 
            std::string name
        )
        {
            std::ofstream outFile;
            outFile.open(
                name.c_str( ), 
                std::ofstream::out | std::ostream::trunc
            );
            if ( !outFile )
            {
                std::cerr << "Can't open file [" << name << "] for output, disable plugin output. " << std::endl;
                isMaster = false; // no Master anymore -> no process is able to write
            }
            else
            {
                for (
                    unsigned int index_direction = 0; 
                    index_direction < parameters::N_observer; 
                    ++index_direction
                ) // over all directions
                {
                    for (
                        unsigned index_omega = 0; 
                        index_omega < radiation_frequencies::N_omega; 
                        ++index_omega
                    ) // over all frequencies
                    {
                        // Take Amplitude for one direction and frequency,
                        // calculate the square of the absolute value
                        // and write to file.
                        constexpr float_X transRadUnit = 
                            UNIT_CHARGE * UNIT_CHARGE * 
                            ( 1.0 / ( 4 * PI * SI::EPS0_SI * PI * PI * SI::SPEED_OF_LIGHT_SI 
                            // * double(particles::TYPICAL_NUM_PARTICLES_PER_MACROPARTICLE) 
                            // * double(particles::TYPICAL_NUM_PARTICLES_PER_MACROPARTICLE) 
                            ) );
                        outFile <<
                            values[
                                index_direction * radiation_frequencies::N_omega + index_omega
                            ] * transRadUnit << "\t";

                    }
                    outFile << std::endl;
                }
                outFile.flush( );
                outFile << std::endl; //now all data are written to file

                if ( outFile.fail( ) )
                    std::cerr << "Error on flushing file [" << name << "]. " << std::endl;

                outFile.close( );
            }
        }

        void 
        calculateRadiationParticles(
            uint32_t currentStep
        )
        {
            DataConnector &dc = Environment< >::get( ).DataConnector( );
            auto particles = dc.get< T_ParticlesType >( 
                T_ParticlesType::FrameType::getName( ), 
                true
            );

            /* execute the particle filter */
            radiation::executeParticleFilter( particles, currentStep );

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
            DataSpace< simDim > localSize( cellDescription->getGridLayout( ).getDataSpaceWithoutGuarding( ) );
            const uint32_t numSlides = MovingWindow::getInstance( ).getSlideCounter( currentStep );
            const SubGrid< simDim >& subGrid = Environment< simDim >::get( ).SubGrid( );
            DataSpace< simDim > globalOffset( subGrid.getLocalDomain( ).offset );
            globalOffset.y( ) += ( localSize.y( ) * numSlides );

            constexpr uint32_t numWorkers = pmacc::traits::GetNumWorkers<
                pmacc::math::CT::volume< SuperCellSize >::type::value
            >::value;

            // PIC-like kernel call of the radiation kernel
            PMACC_KERNEL( KernelTransRadParticles<
                numWorkers
            >{ } )(
                gridDim_rad,
                numWorkers
            )(
                /*Pointer to particles memory on the device*/
                particles->getDeviceParticlesBox( ),

                /*Pointer to memory of radiated amplitude on the device*/
                incTransRad->getDeviceBuffer( ).getDataBox( ),
                cohTransRadPara->getDeviceBuffer( ).getDataBox( ),
                cohTransRadPerp->getDeviceBuffer( ).getDataBox( ),
                numParticles->getDeviceBuffer( ).getDataBox( ),
                globalOffset,
                *cellDescription,
                freqFkt,
                subGrid.getGlobalDomain( ).size
            );

            dc.releaseData( T_ParticlesType::FrameType::getName( ) );
        }
    };

} // namespace picongpu