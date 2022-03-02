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
                fftw_complex *

                // Array to store the result
                vec2r shadowgram;
            public:
                Helper()
                {
                }

                ~Helper()
                {
                }
                
                /** Add the fields from simulation to time integrated fields in Fourier space, is supposed to be called each step the plugin gets called
                 *
                 */
                void append_fields(vec2r Ex, vec2r Ey, vec2r Bx, vec2r By)
                {
                    // Go through all 4 fields
                    fourierEx += propagate_dft(transversal_fft(Ex));
                    fourierEy += propagate_dft(transversal_fft(Ey));
                    fourierBx += propagate_dft(transversal_fft(Bx));
                    fourierBy += propagate_dft(transversal_fft(By));
                }

                /** Calculate the shadowgrams to be called after last simulation step which adds to plugin
                 *
                 */
                void calculate_shadowgrams()
                {
                    // Go through all 4 fields and go back to position space
                    vec3c const tmpEx = longitudinal_idft(transversal_ifft(apply_masks(fourierEx)));
                    vec3c const tmpEy = longitudinal_idft(transversal_ifft(apply_masks(fourierEy)));
                    vec3c const tmpBx = longitudinal_idft(transversal_ifft(apply_masks(fourierBx)));
                    vec3c const tmpBy = longitudinal_idft(transversal_ifft(apply_masks(fourierBy)));

                    perform_time_integration(calculate_poynting_vectors(tmpEx, tmpEy, tmpBx, tmpBy));
                }

                /** Get shadowgram
                 *
                 */
                const vec2r get_shadowgram() const
                {
                    return shadowgram;
                }

            private:

            }
        }
    }
}