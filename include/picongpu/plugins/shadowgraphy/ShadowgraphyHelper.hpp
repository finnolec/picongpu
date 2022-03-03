 #pragma once

#include <fftw3.h>

namespace picongpu
{
    namespace plugins
    {
        namespace shadowgraphy
        {
            class Helper
            {
            private:
                // Arrays for the energy densities E_x * B_y and E_y * B_x
                fftw_complex *;

                // Array to store the result
                vec2r shadowgram;
            public:
                // Constructor of the shadowgraphy helper class
                // To be called at the first time step when the shadowgraphy time integration starts
                Helper()
                {
                    // Create fftw plan for transverse fft for real to complex
                    // Many ffts will be performed -> use FFTW_MEASURE as flag

                    // Create fftw plan for transverse ifft for complex to complex
                    // Even more iffts will be performed -> use FFTW_MEASURE as flag
                }

                // Destructor of the shadowgraphy helper class
                // To be called at the last time step when the shadowgraphy time integration ends
                ~Helper()
                {
                }
                
                // Energy flux calculation loop
                void calculate_energy_flux()
                {
                    // Transversal FFT of E and B fields to get k_x and k_y components
                    // Use fftw plan for fft
                    // E(x, y, zs, tn), B(x, y, zs, tn) -> E(kx, ky, zs, tn), B(kx, ky, zs, tn)

                    // Loop through all omega
                    // E(kx, ky, zs, tn), B(kx, ky, zs, tn) -> E(kx, ky, zs, omega), B(kx, ky, zs, omega)

                        // Apply masks and propagate
                        // E(kx, ky, zs, omega), B(kx, ky, zs, omega) -> M(kx, ky, omega)*E(kx, ky, zo, omega), M(kx, ky omega)*B(kx, ky, zo, omega)

                        // iFFT back into position space
                        // Use fftw plan ifft
                        // M(kx, ky, omega)*E(kx, ky, zo, omega), M(kx, ky, omega)*B(kx, ky, zo, omega) -> E'(x, y, zo, omega), B'(x, y, zo, omega) 

                        // Sum energy fluxes
                        // EF(x, y, zo, omega, tn)
                        // = E'(x, y, zo, omega) * B'(x, y, zo, omega)
                        // + E'(x, y, zo, omega) * Bsum(x, y, zo, omega, tn-1)
                        // + Esum(x, y, zo, omega, tn-1) * B'(x, y, zo, omega)
                        // + EF(x, y, zo, omega, tn-1)

                        // Only do this if it's not the last step of the shadowgraphy integration:
                            // Calculate E edge sum
                            // Esum(x, y, zo, omega, tn) = Esum(x, y, zo, omega, tn-1) + E'(x, y, zo, omega)

                            // Calculate B edge sum
                            // Bsum(x, y, zo, omega, tn) = Bsum(x, y, zo, omega, tn-1) + B'(x, y, zo, omega)
                        // else if this is the last step of the loop:
                            // Free Esum and Bsum from memory ?

                }

                // Calculate the shadowgram after the independent energy fluxes have been calculated for all time steps
                real2darray calculate_shadowgram()
                {
                    // Initialize 2d array for shadowgram

                    // Loop through all omega
                        // shadowgram(x, y) += real(EF1(x, y, zo, omega, tmax)) - real(EF2(x, y, zo, omega, tmax))

                    // return shadowgram
                }

            private:

            }
        }
    }
}