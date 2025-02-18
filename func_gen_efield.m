function laser_E = Efield(t)
%Gaussion envelope , cosine carrier laser, all in a.u.
%   Input:
%   cycle: laser optical cycle
%   simulate_time : simulation time
%   envelope_center : the center of enveloped laser 
%   max_amplitude : laser electric field max amplitude
%   FWHM : laser full width half maximum
%   Output:
%   Efield(t)
    laser_E = input_parameters.max_amplitude * exp(-(t-input_parameters.envelope_center).^2 ./ (input_parameters.FWHM./(2 * (sqrt(2 * log(2))))).^2 ./ 2) .* cos(2*pi/input_parameters.cycle.*t);
end