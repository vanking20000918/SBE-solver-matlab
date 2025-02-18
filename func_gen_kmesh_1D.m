
function kmesh = func_gen_kmesh_1D(lattice_constant, num_points)
% Generate a 1D k-mesh in equal distance from -pi/a to pi/a where a = lattice constant.
% Input: 
% lattice_constant: single value; the lattice constant.
% num_kpoints: single value; the number of k-points from -pi/a to pi/a; equal distance. Better to be odd number so that it can include the Gamma point.  
% Output:
% 1D array; size = num_kpoints.
    kmesh = linspace(-1 * pi /lattice_constant, 1* pi /lattice_constant, num_points);    
end

