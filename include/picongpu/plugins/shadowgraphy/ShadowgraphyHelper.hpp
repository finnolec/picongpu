 #pragma once

#include <fftw3.h>
#include <pmacc/algorithms/math/defines/pi.hpp>
#include "picongpu/simulation_defines.hpp"
#include <cmath> // what
#include "pmacc/assert.hpp"
#include <stdio.h>

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
                using complex_64 = pmacc::math::Complex<float_64>;
                //using complex_64 = pmacc::math::Complex<float_64>;

                typedef std::vector< std::vector< std::vector< complex_64 > > > vec3c;
                typedef std::vector< std::vector< complex_64 > > vec2c;
                typedef std::vector< complex_64 > vec1c;
                typedef std::vector< std::vector< std::vector< float_64 > > > vec3r;
                typedef std::vector< std::vector< float_64 > > vec2r;
                typedef std::vector< float_64 > vec1r;

                // Arrays to store Ex, Ey, Bx and Bz per time step temporarily
                vec2r tmp_Ex, tmp_Ey;
                vec2r tmp_Bx, tmp_By;

                // Arrays for FFTW
                fftw_complex *fftw_in_f_E; // @TODO: Can this be real? Issue is forward / backward FFT
                fftw_complex *fftw_out_f_E;
                fftw_complex *fftw_in_f_B;
                fftw_complex *fftw_out_f_B;
                fftw_complex *fftw_in_b_E;
                fftw_complex *fftw_out_b_E;
                fftw_complex *fftw_in_b_B;
                fftw_complex *fftw_out_b_B;

                fftw_plan plan_forward_E;
                fftw_plan plan_forward_B;
                fftw_plan plan_backward_E;
                fftw_plan plan_backward_B;

                // Arrays for edge sums of 4 field components
                vec3c edge_Ex;
                vec3c edge_Ey;
                vec3c edge_Bx;
                vec3c edge_By;

                // Arrays for the energy densities E_x * B_y and E_y * B_x
                vec3c energydensity_ExBy;
                vec3c energydensity_EyBx;

                // Size of arrays
                int n_x, n_y, n_omegas;

                // Variables for omega calculations @TODO some initializations and bla
                float dt;
                int nt;

                bool isSlidingWindowEnabled;

            public:
                // Constructor of the shadowgraphy helper class
                // To be called at the first time step when the shadowgraphy time integration starts
                Helper(pmacc::math::Size_t<simDim> globalGridSize):
                    n_x(globalGridSize.x() / params::x_res - 2)
                {
                    // Same amount of omegas as ts 
                    // @TODO int division
                    n_omegas = params::omega_n;
                    
                    dt = params::t_res * SI::DELTA_T_SI;
                    nt = params::t_n / params::t_res;

                    std::cout << "initialized with "<< n_x << ", " << n_y << std::endl;

                    // This is currently not allowed to change during plugin run!
                    isSlidingWindowEnabled = MovingWindow::getInstance().isEnabled();

                    if(isSlidingWindowEnabled){
                        // movingWindowCorrection makes the resulting shadowgram smaller if the moving Window is enabled
                        // The resulting loss in the size of the shadowgram comes from the duration of the time integration
                        float const movingWindowCorrection = nt * dt * float_64(SI::SPEED_OF_LIGHT_SI);
                        n_y = (globalGridSize.y() - movingWindowCorrection / SI::CELL_HEIGHT_SI) / params::y_res - 2;
                        PMACC_ASSERT_MSG(n_y > 0, "n_y must be larger than 0, your moving window goes too fast brrrr");
                    } else {
                        n_y = globalGridSize.y() / params::y_res - 2;
                    }

                    // Make sure spacial grid is even
                    n_x = n_x % 2 == 0 ? n_x : n_x - 1;
                    n_y = n_y % 2 == 0 ? n_y : n_y - 1;

                    // Initialization of storage arrays
                    edge_Ex = vec3c(n_x, vec2c(n_y, vec1c(n_omegas)));
                    edge_Ey = vec3c(n_x, vec2c(n_y, vec1c(n_omegas)));
                    edge_Bx = vec3c(n_x, vec2c(n_y, vec1c(n_omegas)));
                    edge_By = vec3c(n_x, vec2c(n_y, vec1c(n_omegas)));
                    energydensity_ExBy = vec3c(n_x, vec2c(n_y, vec1c(n_omegas)));
                    energydensity_EyBx = vec3c(n_x, vec2c(n_y, vec1c(n_omegas)));

                    tmp_Ex = vec2r(n_x, vec1r(n_y));
                    tmp_Ey = vec2r(n_x, vec1r(n_y));
                    tmp_Bx = vec2r(n_x, vec1r(n_y));
                    tmp_By = vec2r(n_x, vec1r(n_y));

                    init_fftw();
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

                // Store fields in helper class with proper resolution
                template<typename F>
                void store_field(int t, pmacc::container::HostBuffer<float3_64, 2>* fieldBuffer1, pmacc::container::HostBuffer<float3_64, 2>* fieldBuffer2)
                {
                    for(int i = 0; i < n_x; ++i){
                        int const grid_i = i * params::x_res;
                        //std::cout << "i:" << i << std::endl;
                        for(int j = 0; j < n_y; ++j){
                            if(isSlidingWindowEnabled){
                                //int const grid_j = j * params::y_res;
                                float const gridPos = float(j * params::y_res) + SI::SPEED_OF_LIGHT_SI * (nt - t - 1) * dt / (SI::CELL_HEIGHT_SI);
                                float const wr = math::fmod(gridPos, 1.0);
                                float const wl = 1.0 - wr;
                                int const grid_j = math::floor(gridPos);

                                if(F::getName() == "E"){
                                    float_64 const Ex0 = (*(fieldBuffer2->origin()(grid_i, grid_j))).x();
                                    float_64 const Ex1 = (*(fieldBuffer2->origin()(grid_i+1, grid_j))).x();
                                    float_64 const Ex2 = (*(fieldBuffer2->origin()(grid_i+2, grid_j))).x();
                                    tmp_Ex[i][j] = ( wl * Ex0 + Ex1 + wr * Ex2 ) / 2.0; 
                                    float_64 const Ey0 = (*(fieldBuffer2->origin()(grid_i, grid_j))).y();
                                    float_64 const Ey1 = (*(fieldBuffer2->origin()(grid_i, grid_j+1))).y();
                                    float_64 const Ey2 = (*(fieldBuffer2->origin()(grid_i, grid_j+2))).y();
                                    tmp_Ey[i][j] = ( wl * Ey0 + Ey1 + wr * Ey2 ) / 2.0;
                                } else {
                                    float_64 const Bx10 = (*(fieldBuffer1->origin()(grid_i, grid_j))).x();
                                    float_64 const Bx11 = (*(fieldBuffer1->origin()(grid_i, grid_j+1))).x();
                                    float_64 const Bx12 = (*(fieldBuffer1->origin()(grid_i, grid_j+2))).x();
                                    float_64 const Bx20 = (*(fieldBuffer2->origin()(grid_i, grid_j))).x();
                                    float_64 const Bx21 = (*(fieldBuffer2->origin()(grid_i, grid_j+1))).x();
                                    float_64 const Bx22 = (*(fieldBuffer2->origin()(grid_i, grid_j+2))).x();
                                    tmp_Bx[i][j] = ( wl * (Bx10 + Bx20) + Bx11 + Bx21 + wr * (Bx12 + Bx22) ) / 4.0;
                                    float_64 const By10 = (*(fieldBuffer1->origin()(grid_i, grid_j))).x();
                                    float_64 const By11 = (*(fieldBuffer1->origin()(grid_i+1, grid_j))).x();
                                    float_64 const By12 = (*(fieldBuffer1->origin()(grid_i+2, grid_j))).x();
                                    float_64 const By20 = (*(fieldBuffer2->origin()(grid_i, grid_j))).x();
                                    float_64 const By21 = (*(fieldBuffer2->origin()(grid_i+1, grid_j))).x();
                                    float_64 const By22 = (*(fieldBuffer2->origin()(grid_i+2, grid_j))).x();
                                    tmp_By[i][j] = ( wl * (By10 + By20) + By11 + By21 + wr * (By12 + By22) ) / 4.0;
                                }
                            } else {
                                int const grid_j = j * params::y_res;

                                // fix yee offset
                                if(F::getName() == "E"){
                                    tmp_Ex[i][j] = ((*(fieldBuffer2->origin()(grid_i, grid_j))).x() + (*(fieldBuffer2->origin()(grid_i+1, grid_j))).x()) / 2.0; 
                                    tmp_Ey[i][j] = ((*(fieldBuffer2->origin()(grid_i, grid_j))).y() + (*(fieldBuffer2->origin()(grid_i, grid_j+1))).y()) / 2.0;
                                } else {
                                    tmp_Bx[i][j] = ((*(fieldBuffer1->origin()(grid_i, grid_j))).x() + (*(fieldBuffer1->origin()(grid_i, grid_j+1))).x()
                                                    + (*(fieldBuffer2->origin()(grid_i, grid_j))).x() + (*(fieldBuffer2->origin()(grid_i, grid_j+1))).x()) / 4.0;
                                    tmp_By[i][j] = ((*(fieldBuffer1->origin()(grid_i, grid_j))).y() + (*(fieldBuffer1->origin()(grid_i+1, grid_j))).y()
                                                    + (*(fieldBuffer2->origin()(grid_i, grid_j))).y() + (*(fieldBuffer2->origin()(grid_i+1, grid_j))).y()) / 4.0;
                                }
                            }
                        }
                    }
                }
                
                // Energy flux calculation loop
                void calculate_energy_flux(int t, bool is_exby)
                /**
                t: current time step, from 0 to (nt-1)
                is_exby: bool, true if first part of poynting vector, false if second part of poynting vector
                **/
                {
                    // Transversal FFT of E and B fields to get k_x and k_y components
                    // Use fftw plan for fft
                    // E(x, y, zs, tn), B(x, y, zs, tn) -> E(kx, ky, zs, tn), B(kx, ky, zs, tn)
                    for(int i = 0; i < n_x; i++){
                        for(int j = 0; j < n_y; j++){
                            int index = i + j * n_x;

                            // Real values
                            if(is_exby){
                                fftw_in_f_E[index][0] = tmp_Ex[i][j];
                                fftw_in_f_B[index][0] = tmp_By[i][j];
                            } else {
                                fftw_in_f_E[index][0] = tmp_Ey[i][j];
                                fftw_in_f_B[index][0] = tmp_Bx[i][j];
                            }

                            // Imaginary values
                            fftw_in_f_E[index][1] = 0.0;
                            fftw_in_f_B[index][1] = 0.0;
                        }
                    }

                    fftw_execute(plan_forward_E);
                    fftw_execute(plan_forward_B);

                    vec2c E_k(n_x, vec1c(n_y));
                    vec2c B_k(n_x, vec1c(n_y));

                    for(int i = 0; i < n_x; i++){
                        // Put origin into center of array with this, necessary due to FFT
                        int const i_ffs = (i + n_x/2) % n_x;

                        for(int j = 0; j < n_y; j++){
                            // Put origin into center of array with this, necessary due to FFT
                            int const j_ffs = (j + n_y / 2) % n_y;

                            int const index = i + j * n_x;
                            // @TODO write this in nice
                            E_k[i_ffs][j_ffs] = complex_64(fftw_out_f_E[index][0], fftw_out_f_E[index][1]);
                            B_k[i_ffs][j_ffs] = complex_64(fftw_out_f_B[index][0], fftw_out_f_B[index][1]);
                        }
                    }

                    // Loop through all omega
                    // E(kx, ky, zs, tn), B(kx, ky, zs, tn) -> E(kx, ky, zs, omega), B(kx, ky, zs, omega)
                    for(int o = 0; o < n_omegas; o++){
                        // Omega for time domain Fourier trafo
                        //float_64 omega = 2.0 * pmacc::math::Pi<float_64>::value  * (o - nt / 2.0) / nt / dt;
                        float_64 const omega_SI = fourierhelper::omega(o);

                        //printf("185 %e - %e - %e - %e \n", pmacc::math::Pi<float_64>::value, o, nt, dt);

                        // Apply masks and propagate
                        // E(kx, ky, zs, omega), B(kx, ky, zs, omega) -> M(kx, ky, omega)*E(kx, ky, zo, omega), M(kx, ky omega)*B(kx, ky, zo, omega)
                        for(int i = 0; i < n_x; i++){
                            for(int j = 0; j < n_y; j++){
                                int const index = i + j * n_x;

                                float_64 const t_SI = t * int(params::t_res) * float_64(picongpu::SI::DELTA_T_SI);
                                //float_64 const propagator = float_64(params::delta_z) / float_64(SPEED_OF_LIGHT);
                                float_64 const sqrt1 = 1 / (float_64(SI::SPEED_OF_LIGHT_SI) * float_64(SI::SPEED_OF_LIGHT_SI));
                                //printf("c %e \n", SPEED_OF_LIGHT_SI)
                                //printf("1 %e \n", sqrt1);
                                float_64 const sqrt2 = fourierhelper::kx(i) * fourierhelper::kx(j) / (omega_SI * omega_SI);
                                //printf("2 %e \n", sqrt2);
                                float_64 const sqrt3 = fourierhelper::ky(j) * fourierhelper::ky(j) / (omega_SI * omega_SI);
                                //printf("3 %e \n", sqrt3);
                                
                                // für 3 und 4
                                //float_64 const propagator = float_64(params::delta_z) * (picongpu::math::sqrt(sqrt1 - sqrt2 - sqrt3) - omega_SI / float_64(picongpu::SI::DELTA_T_SI));
                                // für 5 und 6 und 7

                                // sqrtContent equal to kz^2 / omega^2
                                float_64 const sqrtContent = sqrt1 - sqrt2 - sqrt3;

                                if (sqrtContent >= 0.0){
                                    float_64 const propagator = float_64(params::delta_z) * (0.0 * picongpu::math::sqrt(sqrtContent) - 1.0 / float_64(SI::SPEED_OF_LIGHT_SI));
                                    //float_64 const propagator = 0 * float_64(params::delta_z) / SI::SPEED_OF_LIGHT_SI;

                                    complex_64 const phase_e = complex_64(0, +omega_SI * (propagator - t_SI));
                                    complex_64 const tmp_e = complex_64(masks::mask(i, j, o)) * E_k[i][j] * math::exp(phase_e);

                                    complex_64 const phase_b = complex_64(0, +omega_SI * (-propagator + t_SI));
                                    complex_64 const tmp_b = complex_64(masks::mask(i, j, o)) * B_k[i][j] * math::exp(phase_b);

                                    fftw_in_b_E[index][0] = tmp_e.get_real();
                                    fftw_in_b_E[index][1] = tmp_e.get_imag();
                                    fftw_in_b_B[index][0] = tmp_b.get_real();
                                    fftw_in_b_B[index][1] = tmp_b.get_imag();
                                } else {
                                    fftw_in_b_E[index][0] = 0.0;
                                    fftw_in_b_E[index][1] = 0.0;
                                    fftw_in_b_B[index][0] = 0.0;
                                    fftw_in_b_B[index][1] = 0.0;
                                }
                            }
                        }

                        // iFFT back into position space
                        // Use fftw plan ifft
                        // M(kx, ky, omega)*E(kx, ky, zo, omega), M(kx, ky, omega)*B(kx, ky, zo, omega) -> E'(x, y, zo, omega), B'(x, y, zo, omega) 
                        fftw_execute(plan_backward_E);
                        fftw_execute(plan_backward_B);

                        for(int i = 0; i < n_x; i++){
                            for(int j = 0; j < n_y; j++){
                                int index = i + j * n_x;

                                complex_64 const E = complex_64(fftw_out_b_E[index][0], fftw_out_b_E[index][1]);
                                complex_64 const B = complex_64(fftw_out_b_B[index][0], fftw_out_b_B[index][1]);


                                if (is_exby)
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
                                        //std::cout << "232 " << edge_Ex[i][j][o].get_real() << ", "<< edge_Ex[i][j][o].get_imag() << std::endl;

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
                                    //std::cout << "248 " << energydensity_EyBx[i][j][o].get_real() << ", "<< energydensity_EyBx[i][j][o].get_imag() << std::endl;

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
                vec2r get_shadowgram() // @TODO: make this return a pointer
                {
                    vec2r shadowgram(n_x, vec1r(n_y));

                    // Loop through all omega
                    for(int i = 0; i < n_x; i++){
                        for(int j = 0; j < n_y; j++){
                            for(int o = 0; o < n_omegas; ++o)
                            {
                                // shadowgram(x, y) += real(EF1(x, y, zo, omega, tmax)) - real(EF2(x, y, zo, omega, tmax))
                                shadowgram[i][j] += (energydensity_ExBy[i][j][o] - energydensity_EyBx[i][j][o]).get_real() * fourierhelper::get_omega_step() / (nt * nt);
                            }
                        }
                    }

                    return shadowgram;
                }

                int get_n_x(){
                    return n_x;
                }

                int get_n_y(){
                    return n_y;
                }

            private:
                // Initialize fftw memory things, supposed to be called once per plugin loop
                void init_fftw()
                {
                    std::cout << "init fftw" << std::endl;
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
                    plan_forward_E = fftw_plan_dft_2d(n_y, n_x, fftw_in_f_E, fftw_out_f_E, FFTW_FORWARD, FFTW_MEASURE);
                    plan_forward_B = fftw_plan_dft_2d(n_y, n_x, fftw_in_f_B, fftw_out_f_B, FFTW_BACKWARD, FFTW_MEASURE);
                    //plan_forward_B = fftw_plan_dft_2d(n_y, n_x, fftw_in_f_B, fftw_out_f_B, FFTW_FORWARD, FFTW_MEASURE);

                    // Create fftw plan for transverse ifft for complex to complex
                    // Even more iffts will be performed -> use FFTW_MEASURE as flag
                    plan_backward_E = fftw_plan_dft_2d(n_y, n_x, fftw_in_b_E, fftw_out_b_E, FFTW_BACKWARD, FFTW_MEASURE);
                    plan_backward_B = fftw_plan_dft_2d(n_y, n_x, fftw_in_b_B, fftw_out_b_B, FFTW_FORWARD, FFTW_MEASURE);
                    //plan_backward_B = fftw_plan_dft_2d(n_y, n_x, fftw_in_b_B, fftw_out_b_B, FFTW_BACKWARD, FFTW_MEASURE);
                }

            }; // class Helper
        } // namespace shadowgraphy
    } // namespace plugins
} //namespace picongpu