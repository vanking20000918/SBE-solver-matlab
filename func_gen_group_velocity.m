function velocity = func_gen_group_velocity(kmesh, band)
%Calculate group velocity according k grids and bandstructure.
%   Input: k grids, bandstructure
%   Output: group velocity, a.u.
dk = kmesh(2) - kmesh(1);
velocity = zeros(1, length(kmesh));
velocity(1) = (band(2) - band(end-1)) / 2/ dk;
for i = 2:(length(kmesh)-1)
    velocity(i) = (band(i+1)-band(i-1))/2/dk;
end
velocity(end) = velocity(1);
end