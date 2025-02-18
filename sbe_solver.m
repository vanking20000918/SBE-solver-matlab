% 1D 2 bands SBEs
clearvars -except files;clc;tic % time begins
disp('Solver begins')
%% Part 1: tight-banding model  
% Bandstructure & group velocity
kmesh = func_gen_kmesh_1D(input_parameters.lattice_constant, input_parameters.num_kpoints);
for kpoint_index = 1:length(kmesh)  
    conduction_band(kpoint_index) = func_gen_tight_binding_band_1D(input_parameters.lattice_constant, kmesh(kpoint_index), input_parameters.conduction_band_param, input_parameters.conduction_band_shift); % gen CB  
    valence_band(kpoint_index) = func_gen_tight_binding_band_1D(input_parameters.lattice_constant, kmesh(kpoint_index), input_parameters.valence_band_param, input_parameters.valence_band_shift); % gen VB  
end  
conduction_velocity = func_gen_group_velocity(kmesh, conduction_band);% Periodic boundary condition
valence_velocity = func_gen_group_velocity(kmesh, valence_band);
% plot bandstructure
figure(1)
plot(kmesh*input_parameters.lattice_constant/2/pi, conduction_band* 27.211386245988, color = 'r');
hold on
plot(kmesh*input_parameters.lattice_constant/2/pi, valence_band*27.211386245988, color = 'b');
hold off
title('Tight-binding bands')
xlabel('k (reciprocal vector)');
ylabel('Energy (eV)');
legend('CB', 'VB');
exportgraphics(gcf,'Tight-binding bands.png','Resolution',300)
% plot group vilocity
figure(2)
plot(kmesh*input_parameters.lattice_constant/2/pi, conduction_velocity, color = 'r');
hold on
plot(kmesh*input_parameters.lattice_constant/2/pi, valence_velocity, color = 'b');
hold off
title('Band group velocities')
xlabel('k (reciprocal vector)');
ylabel('Group velocity (a.u.)');
legend('CB', 'VB');
exportgraphics(gcf,'Band group velocities.png','Resolution',300);

%% Part 2 laser
tmesh = linspace(0, input_parameters.simulate_time, input_parameters.num_time_step);
laser_E = func_gen_efield(tmesh);
for t_index = 1:length(tmesh)
    laser_A(t_index) = -1 * integral(@func_gen_efield, tmesh(1), tmesh(t_index));
end
% plot Efield
figure(3)
plot(tmesh/ input_parameters.cycle, laser_E);
title('Laser E field')
xlabel('Time (cycles)');
ylabel('Efield (a.u.)');
exportgraphics(gcf,'Laser E field.png','Resolution',300);
hold off
% plot Afield
figure(4)
plot(tmesh/ input_parameters.cycle, laser_A);
title('Laser A field')
xlabel('Time (cycles)');
ylabel('Afield (a.u.)');
exportgraphics(gcf,'Laser A field.png','Resolution',300);
hold off
%% Part 3: k-dependent transition dipole moment
if exist("transition dipole moment.txt", "file") == 2
    disp("File 'transition dipole moment.txt' exist.")
    if length(readmatrix("transition dipole moment.txt")) ~= length(kmesh)
        error("Size does not match!")
        return
    else
        TDM = func_gen_transition_dipole_moment_1D(kmesh); % k-dependent transition dipole moment
        disp("Read it successfully!")
    end
else
    error("File 'transition dipole moment' does not exist.")
    return
end

if length(TDM) ~= length(kmesh)
    disp("Size doesn't match")
end
% plot transition dipole moment
figure(5)
plot(kmesh*input_parameters.lattice_constant/2/pi, TDM, color = 'r');
title('Transition dipole moment');
xlabel('k (reciprocal vector)');
ylabel('Transition dipole moment(a.u.)');
exportgraphics(gcf, 'Transition dipole moment.png','Resolution',300);
hold off
%% Part 4: numerical solution to SBEs
% Warning! This part is very time-consuming!!!
x=kmesh;
t=tmesh;
m=0;
options = odeset('RelTol',1e-6,'AbsTol',1e-6,'InitialStep',10);
if exist("solution.mat", "file") == 2
    disp("'solution.mat' exist!")
    read_variables = input("Load it:Y/N?   ","s");
     if (read_variables == "Y") || (read_variables == "y")
        load("solution.mat");
     else
         sol=pdepe(m,@pdefun,@pdeic,@pdebc,x,t,options);
     end
else 
    sol=pdepe(m,@pdefun,@pdeic,@pdebc,x,t,options);
end

p=sol(:,:,1)/input_parameters.num_kpoints;
fe=real(sol(:,:,2))/input_parameters.num_kpoints;
fh=real(sol(:,:,3))/input_parameters.num_kpoints;
p_tol = sum(p,2); % micropolarization
fe_tol = sum(fe,2); % occupation of electrons in CB
fh_tol = sum(fh,2); % occupation of holes in VB
% plot p, fe, fh
figure(6)
subplot(2,1,1);
plot(t/input_parameters.cycle,real(p_tol))
title('Real part of micropolarization')
xlabel('Time(cycles)')
ylabel('Real(p)')
subplot(2,1,2);
plot(t/input_parameters.cycle,imag(p_tol))
title('Imaginary part of micropolarization')
xlabel('Time(cycles)')
ylabel('Imag(p)')
exportgraphics(gcf, 'Micropolarization.png','Resolution',300);
hold off
figure(7)
subplot(2,1,1)
plot(t/input_parameters.cycle,fe_tol, color = 'red')
title('Electrons occupation of CB')
xlabel('Time(cycles)')
ylabel('fe(relative value)')
subplot(2,1,2)
plot(t/input_parameters.cycle,fh_tol, color = 'blue')
title('Holes occupation of VB')
xlabel('Time(cycles)')
ylabel('fh(relative value)')
exportgraphics(gcf, 'fe & fh.png','Resolution',300);
hold off
%% Part 5: calculate the currents & HHG
P = sum(real(TDM' .* p) *2, 2)'; % macroscropic polarization, note big 'P'
Jer = zeros(1, length(tmesh));
for t_index = 2:(length(tmesh)-1)
    Jer(t_index) = (P(t_index+1) - P(t_index-1))/2/(tmesh(2)-tmesh(1));
end
Jra = sum(conduction_velocity .* fe + (-1) * valence_velocity .* fh, 2)'; % intraband current
Jtot = Jer + Jra; % total current
% plot interband & intraband currents
figure(8)
subplot(3,1,1)
plot(tmesh/ input_parameters.cycle, Jer, color = 'blue')
title('Interband current')
xlabel('Time(cycles)')
ylabel('Jer(a.u.)')
subplot(3,1,2)
plot(tmesh/ input_parameters.cycle, Jra, color = 'red')
title('Intraband current')
xlabel('Time(cycles)')
ylabel('Jra(a.u.)')
subplot(3,1,3)
plot(tmesh/ input_parameters.cycle, Jtot, color = 'black')
title('Total current')
xlabel('Time(cycles)')
ylabel('Jtot(a.u.)')
exportgraphics(gcf, 'Currents.png','Resolution',300);
hold off
% FFT
h = 6.62607015e-34; % Planck constant, J*s 
e = 1.602176634e-19; % elementary charge, C
dt = input_parameters.simulate_time / (input_parameters.num_time_step - 1); % time interval, a.u.
L = length(tmesh);
window = hann(L)'; % available window function: hamming, bartlett, blackman, gausswin, hann, so on.
Jer_window = Jer .* window; % window processing
Jra_window = Jra .* window;
Jtot_window = Jtot .* window;
if mod(length(tmesh),2) ~= 0 % convert length(tmesh) to even
    L = L -1;
end
dfrqes = 1/dt; % frquency interval, a.u.

coefficient_er = fft(Jer_window); % FFT coefficient
coefficient_ra = fft(Jra_window);
coefficient_tot = fft(Jtot_window);
freqs = dfrqes * (0:(L/2))/L;  % generate shifted frequency sequence, a.u.
freqs = freqs .* 41.341374575751 .* (h/1e-15/e) / input_parameters.photon_energy; % unit convertion from a.u. to photon energy
double_spectrum_er = abs(coefficient_er/L);
sigle_spectrum_er = double_spectrum_er(1:L/2+1);
sigle_spectrum_er(2:end-1) = 2* sigle_spectrum_er(2:end-1);
double_spectrum_ra = abs(coefficient_ra/L);
sigle_spectrum_ra = double_spectrum_ra(1:L/2+1);
sigle_spectrum_ra(2:end-1) = 2* sigle_spectrum_ra(2:end-1);
double_spectrum_tot = abs(coefficient_tot/L);
sigle_spectrum_tot = double_spectrum_tot(1:L/2+1);
sigle_spectrum_tot(2:end-1) = 2* sigle_spectrum_tot(2:end-1);

Power_er = abs(sigle_spectrum_er) .^2; % FFT power
Power_ra = abs(sigle_spectrum_ra) .^2;
Power_tot = abs(sigle_spectrum_tot) .^2;
% plot HHG of jer, jra, jtot
figure(9)
subplot(3,1,1)
plot(freqs, log(Power_er), color = 'blue')
xlim([1, input_parameters.Order+1])
ylim([-inf, inf])
set(gca, 'XTick', 1:10:(input_parameters.Order+1)) % set xtick reference line
title('Interband HHG')
xlabel('Frequency(Photon energy)')
ylabel('Log(power)(arb.units)')
subplot(3,1,2)
plot(freqs, log(Power_ra), color = 'red')
xlim([1, input_parameters.Order+1])
ylim([-inf, inf])
set(gca, 'XTick', 1:10:(input_parameters.Order+1))
title('Intraband HHG')
xlabel('Frequency(Photon energy)')
ylabel('Log(power)(arb.units)')
subplot(3,1,3)
plot(freqs, log(Power_tot), color = 'black')
xlim([1, input_parameters.Order+1])
ylim([-inf, inf])
set(gca, 'XTick', 1:10:(input_parameters.Order+1))
title('Total HHG')
xlabel('Frequency(Photon energy)')
ylabel('Log(power)(arb.units)')
exportgraphics(gcf, 'Fourier transform analysis.png','Resolution',300);
hold off

save('solution.mat', "sol") % save variable 'sol'
save("All variables.mat") % save all variables 
disp('Solver completed')
toc % time up
%% Part 6: semiconductor bloch equations
function [c, f, s]=pdefun(k, t, u, uk) % u(1), u(2), u(3) representing micropolarization, occupation of electrons in CB, occupation of holes in VB.
    c=[1i; 1; 1] ;
    conduction_band = func_gen_tight_binding_band_1D(input_parameters.lattice_constant, k, input_parameters.conduction_band_param, input_parameters.conduction_band_shift); % CB
    valence_band = func_gen_tight_binding_band_1D(input_parameters.lattice_constant, k, input_parameters.valence_band_param, input_parameters.valence_band_shift); % VB
    epsilon_e = conduction_band; % omitting coulomb interaction
    epsilon_h = -1 * valence_band;
    T2 = input_parameters.T2; % dephasing time of polarization
    transition_dipole_moment = func_gen_transition_dipole_moment_1D(k); % transition dipole moment
    laser_E = func_gen_efield(t); % Efield 
    s1 = (epsilon_h+epsilon_e-1i /T2)*u(1) - (1-u(2)- u(3))*transition_dipole_moment*laser_E; % s1-coefficient
    s2 = -2* imag(transition_dipole_moment* laser_E * conj(u(1))); % s2-coefficient
    s3 =  s2; % s3-coefficient
    s= [s1; s2; s3];
    f= laser_E * [1i*u(1); u(2); u(3)];
    
    timebar = t / input_parameters.simulate_time; % show waitbar
    percent_timebar = strcat(num2str(timebar*100),'%');
    
    fprintf('Solving, %s.\n', percent_timebar)
end

%Initial condition
function [u0]=pdeic(k)
    u0=[0; 0; 0];
end

%Boundary condition
% ps: we fix the micropolarization and occupation of boundary k points, so
% you need enough numbers of kpoints(num_points).
% k = -pi/a and pi/a, p = 0, fe= 0, fh = 0. 
function [pl, ql, pr, qr]=pdebc(xl ,ul ,xr ,ur, t)
    %laser_E = func_gen_efield(t); % Efield 
    pl= ul;
    ql=[0; 0; 0];
    pr= ur;
    qr=[0; 0; 0];
end
