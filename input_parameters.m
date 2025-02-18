classdef input_parameters < handle

    properties(Constant)

       % Input parameters, all units are in a.u.
       num_kpoints = 51; % number of k points; better to be odd number so that it can include the Gamma point.
       
       lattice_constant = 5.32; % lattice constant;
       conduction_band_param = [0.0898, -0.0814, -0.0024, -0.0048, -0.0003, -0.0009]; % CB expansion coefficients
       valence_band_param = [-0.0928, 0.0705, 0.0200, -0.0012, 0.0029, 0.0006]; % VB expansion coefficients
       conduction_band_shift = 3.3 / 27.211386245988; % CB Energy shift
       valence_band_shift = 0 / 27.211386245988; % VB Energy shift
       
       cycle = 10.9 * 41.341374575751; % laser optical cycle
       photon_energy = 2* pi/ input_parameters.cycle * 27.211386245988; % laser photon energy, eV
       simulate_time = 5 * input_parameters.cycle; % simulation time
       num_time_step = 546; % dt = simulate_time / (num_time_step - 1)
       envelope_center = input_parameters.simulate_time / 2; % the center of enveloped laser 
       max_amplitude = 0.003; % laser electric field max amplitude
       FWHM = 2 * input_parameters.cycle; % laser full width half maximum
       T2 = input_parameters.cycle /4; % polarization's dephasing time
       %T1 = input_parameters.simulate_time * 100 % carriers occupation's dephasing time  
       Order = 40; % order of FFT
    end
end

