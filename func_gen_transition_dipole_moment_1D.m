function TDM = func_gen_transition_dipole_moment_1D(kpoint_coord)
%   Fit transition dipole moment corresponding to k point_coord by parabola.
%   Input:
%   kpoint_coord: single value; the reciprocal coordinate of kpoint.
%   Output:
%   transition dipole moment, a.u.
transition_dipole_moment = readmatrix("transition dipole moment.txt"); % transition dipole moment,a.u., read from file "transition dipole moment.txt"
TDM_fit = fit(transpose(func_gen_kmesh_1D(input_parameters.lattice_constant, input_parameters.num_kpoints)), transition_dipole_moment, 'poly2');
TDM = TDM_fit(kpoint_coord);
end