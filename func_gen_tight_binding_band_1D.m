function energy = func_gen_tight_binding_band_1D(lattice_constant, kpoint_coord, gen_param, energy_shift)
% Generate band energy at given kpoint from cosine function expansion.
% Input:
% lattice_constant: single value; the lattice constant.
% kpoint_coord: single value; the reciprocal coordinate of kpoint.
% gen_param: 1D numpy array; the j-th coefficient of cos( j * k * a ) where a = lattice constant.
% energy_shift: single value; extra energy added on the cosine expansion; e.g., the band gap.
% Output:
% The band energy as a single value
    energy = 0;  
    for param_index = 1:length(gen_param)  
        energy = energy + gen_param(param_index) * cos((param_index-1) * kpoint_coord * lattice_constant);  
    end
    energy = energy + energy_shift;
end