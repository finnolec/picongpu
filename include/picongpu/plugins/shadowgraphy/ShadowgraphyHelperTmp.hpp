 #pragma once

#include <fftw3.h>

namespace picongpu
{
    namespace plugins
    {
        namespace shadowgraphy
        {
            using namespace std::complex_literals;

            class Helper
            {
            private:
                typedef std::vector< std::vector< std::vector< complex_64 > > > vec3c;
                typedef std::vector< std::vector< complex_64 > > vec2c;
                typedef std::vector< complex_64 > vec1c;
                typedef std::vector< std::vector< std::vector< float_X > > > vec3r;
                typedef std::vector< std::vector< float_X > > vec2r;

                // Arrays for FFTW
                fftw_complex *fftw_in_f_E; // @TODO: Can this be real? Issue is forward / backward FFT
                fftw_complex *fftw_out_f_E;
                fftw_complex *fftw_in_f_B;
                fftw_complex *fftw_out_f_B;
                fftw_complex *fftw_in_b_E;
                fftw_complex *fftw_out_b_E;
                fftw_complex *fftw_in_b_B;
                fftw_complex *fftw_out_b_B;

                // Arrays for edge sums of 4 field components
                vec3c64 edge_Ex;
                vec3c64 edge_Ey;
                vec3c64 edge_Bx;
                vec3c64 edge_By;

                // Arrays for the energy densities E_x * B_y and E_y * B_x
                vec3c64 energydensity_ExBy;
                vec3c64 energydenstiy_EyBx;

                // Size of arrays
                int n_x, n_y, n_omegas;

                // Variables for omega calculations @TODO some initializations and bla
                float dt;
                int nt;

                float delta_z;

            public:
                // Constructor of the shadowgraphy helper class
                // To be called at the first time step when the shadowgraphy time integration starts
                Helper(int n_x, int n_y, float delta_z):
                    n_x(n_x),
                    n_y(n_y),
                    delta_z(delta_z)
                {
                    // Input and output arrays for the FFT transforms
                    fftw_in_f_E = fftw_alloc_complex(n_x * n_y);
                    fftw_out_f_E = fftw_alloc_complex(n_x * n_y);
                    fftw_in_f_B = fftw_alloc_complex(n_x * n_y);
                    fftw_out_f_B = fftw_alloc_complex(n_x * n_y);

                    fftw_in_b_E = fftw_alloc_complex(n_x * n_y);
                    fftw_out_b_E = fftw_alloc_complex(n_x * n_y);
                    fftw_in_b_B = fftw_alloc_complex(n_x * n_y);
                    fftw_out_b_B = fftw_alloc_complex(n_x * n_y);

                    // Create fftw plan for transverse fft for real to complex
                    // Many ffts will be performed -> use FFTW_MEASURE as flag
                    plan_forward_E = fftw_plan_dft_2d(n_x, n_y, fftw_in_f_E, fftw_out_f_E, FFTW_FORWARD, FFTW_MEASURE);
                    plan_forward_B = fftw_plan_dft_2d(n_x, n_y, fftw_in_f_B, fftw_out_f_B, FFTW_BACKWARD, FFTW_MEASURE);

                    // Create fftw plan for transverse ifft for complex to complex
                    // Even more iffts will be performed -> use FFTW_MEASURE as flag
                    plan_backward_E = fftw_plan_dft_2d(n_x, n_y, fftw_in_b_E, fftw_out_b_E, FFTW_BACKWARD, FFTW_MEASURE);
                    plan_packward_B = fftw_plan_dft_2d(n_x, n_y, fftw_in_b_B, fftw_out_b, FFTW_FORWARD, FFTW_MEASURE);

                    // Initialization of storage arrays
                    edge_Ex = vec3c(n_x, vec2c(n,y, vec1c(n_omegas)));
                    edge_Ey = vec3c(n_x, vec2c(n,y, vec1c(n_omegas)));
                    edge_Bx = vec3c(n_x, vec2c(n,y, vec1c(n_omegas)));
                    edge_By = vec3c(n_x, vec2c(n,y, vec1c(n_omegas)));
                    energydensity_ExBy = vec3c(n_x, vec2c(n,y, vec1c(n_omegas)));
                    energydenstiy_EyBx = vec3c(n_x, vec2c(n,y, vec1c(n_omegas)));
                }

                // Destructor of the shadowgraphy helper class
                // To be called at the last time step when the shadowgraphy time integration ends
                ~Helper()
                {
                    fftw_free(fftw_in_f_E);
                    fftw_free(fftw_out_f_E);
                    fftw_free(fftw_in_f_B);
                    fftw_free(fftw_out_f_B);
                    fftw_free(fftw_in_b_E);
                    fftw_free(fftw_out_b_E);
                    fftw_free(fftw_in_b_B);
                    fftw_free(fftw_out_b_B);
                }
                
                // Energy flux calculation loop
                template<typename T>
                void calculate_energy_flux(T const& sim_E, T const& sim_B, int t, bool is_first_summand)
                /**
                sim_E: E from simulation, real 2d array
                sim_B: B from simulation, real 2d array
                t: current time step, from 0 to (nt-1)
                is_first_summand: bool, true if first part of poynting vector, false if second part of poynting vector
                **/
                {
                    
                    // Transversal FFT of E and B fields to get k_x and k_y components
                    // Use fftw plan for fft
                    // E(x, y, zs, tn), B(x, y, zs, tn) -> E(kx, ky, zs, tn), B(kx, ky, zs, tn)
                    for(int i = 0; i < n_x; i++){
                        for(int j = 0; j < n_y; j++){
                            fftw_in_f_E[i][j][0] = sim_E[i][j];
                            fftw_in_f_E[i][j][1] = 0.0;

                            fftw_in_f_B[i][j][0] = sim_B[i][j];
                            fftw_in_f_B[i][j][1] = 0.0;
                        }
                    }

                    fftw_execute(plan_forward_E);
                    fftw_execute(plan_forward_B);

                    vec2c E_k(n_x, vec1c(n_y));
                    vec2c B_k(n_x, vec1c(n_y));

                    for(int i = 0; i < n_x; i++){
                        for(int j = 0; j < n_y; j++){
                            // @TODO write this in nice
                            E_k[i][j] = fftw_out_f_E[i][j][0] + 1i * fftw_out_f_E[i][j][1];

                            B_k[i][j] = fftw_out_f_B[i][j][0] + 1i * fftw_out_f_B[i][j][1];
                        }
                    }

                    // Loop through all omega
                    // E(kx, ky, zs, tn), B(kx, ky, zs, tn) -> E(kx, ky, zs, omega), B(kx, ky, zs, omega)
                    for(int o = 0; o < n_omegas; o++){
                        // Omega for time domain Fourier trafo
                        float_X omega = 2 * math::PI  * (o - nt / 2.0) / nt / dt;

                        // Apply masks and propagate
                        // E(kx, ky, zs, omega), B(kx, ky, zs, omega) -> M(kx, ky, omega)*E(kx, ky, zo, omega), M(kx, ky omega)*B(kx, ky, zo, omega)
                        for(int i = 0; i < n_x; i++){
                            for(int j = 0; j < n_y; j++){
                                fftw_in_b_E[i][j] = mask(i, j, o) * E_k[i][j] * math::exp(-1i * omega * (t - delta_z / const::c));
                                fftw_in_b_B[i][j] = mask(i, j, o) * B_k[i][j] * math::exp(+1i * omega * (t + delta_z / const::c));
                            }
                        }

                        // iFFT back into position space
                        // Use fftw plan ifft
                        // M(kx, ky, omega)*E(kx, ky, zo, omega), M(kx, ky, omega)*B(kx, ky, zo, omega) -> E'(x, y, zo, omega), B'(x, y, zo, omega) 
                        fftw_execute(plan_backward_E);
                        fftw_execute(plan_backward_B);

                        for(int i = 0; i < n_x; i++){
                            for(int j = 0; j < n_y; j++){
                                constexpr complex_X E = fftw_out_b_E[i][j][0] + 1i * fftw_out_b_E[i][j][1];
                                constexpr complex_X B = fftw_out_b_B[i][j][0] + 1i * fftw_out_b_B[i][j][1];

                                if (is_first_summand)
                                {
                                    // Sum energy fluxes
                                    // EF(x, y, zo, omega, tn)
                                    // = E'(x, y, zo, omega) * B'(x, y, zo, omega)
                                    // + E'(x, y, zo, omega) * Bsum(x, y, zo, omega, tn-1)
                                    // + Esum(x, y, zo, omega, tn-1) * B'(x, y, zo, omega)
                                    // + EF(x, y, zo, omega, tn-1)
                                    energydensity_ExBy[i][j][o] += E * B + edge_Ex[i][j][o] * B + E * edge_By[i][j][o];

                                    // Only do this if it's not the last step of the shadowgraphy integration:
                                    if (t < ( nt - 1 )) {
                                        // Calculate E edge sum
                                        // Esum(x, y, zo, omega, tn) = Esum(x, y, zo, omega, tn-1) + E'(x, y, zo, omega)
                                        edge_Ex[i][j][o] += E;

                                        // Calculate B edge sum
                                        // Bsum(x, y, zo, omega, tn) = Bsum(x, y, zo, omega, tn-1) + B'(x, y, zo, omega)
                                        edge_By[i][j][o] += B;
                                    }
                                    // else if this is the last step of the loop: @TODO
                                        // Free Esum and Bsum from memory ?
                                } else {
                                    // Sum energy fluxes
                                    // EF(x, y, zo, omega, tn)
                                    // = E'(x, y, zo, omega) * B'(x, y, zo, omega)
                                    // + E'(x, y, zo, omega) * Bsum(x, y, zo, omega, tn-1)
                                    // + Esum(x, y, zo, omega, tn-1) * B'(x, y, zo, omega)
                                    // + EF(x, y, zo, omega, tn-1)
                                    energydensity_EyBx[i][j][o] += E * B + edge_Ey[i][j][o] * B + E * edge_Bx[i][j][o];

                                    // Only do this if it's not the last step of the shadowgraphy integration: 
                                    if (t < ( nt - 1 )) {
                                        // Calculate E edge sum
                                        // Esum(x, y, zo, omega, tn) = Esum(x, y, zo, omega, tn-1) + E'(x, y, zo, omega)
                                        edge_Ey[i][j][o] += E;

                                        // Calculate B edge sum
                                        // Bsum(x, y, zo, omega, tn) = Bsum(x, y, zo, omega, tn-1) + B'(x, y, zo, omega)
                                        edge_Bx[i][j][o] += B;
                                    }
                                    // else if this is the last step of the loop: @TODO
                                        // Free Esum and Bsum from memory ?

                                }
                            }
                        }

                        

                    }
                }

                // Calculate the shadowgram after the independent energy fluxes have been calculated for all time steps
                vec2r calculate_shadowgram()
                {
                    vec2r shadowgram(n_x, vec1c(n_y));

                    // Loop through all omega
                    for(int i = 0; i < n_x; i++){
                        for(int j = 0; j < n_y; j++){
                            for(int o = 0; o < n_omegas; ++o)
                            {
                                // shadowgram(x, y) += real(EF1(x, y, zo, omega, tmax)) - real(EF2(x, y, zo, omega, tmax))
                                shadowgram[i][j] += (energydensity_ExBy[i][j][o] - energydensity_EyBx[i][j][o]).real();
                            }
                        }
                    }

                    return shadowgram;
                }

            private:

            }
        }
    }
}