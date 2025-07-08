% -----------
% Case 01
% -----------
% To make life easy, the surface area is 1 m2.
% The thickness is set by parameter dz.
% Thus the volume is 1 x dz m3.
%------------
%
%--------------------------------------------------------------------
% In case 01.1 the box is in equilibrium with the atmosphere before 
% the alkalinity perturbation at t=0.
%--------------------------------------------------------------------

% All relevant parameters are defined in a structure called case_**
% that can be passed on to the functions that do the calculations.

clear; close all;

plot_output = 1; % yes/no

% INPUT PARAMETERS 

% time-stepping parameters
secpday   = 24*3600;        % seconds per day 
Tmax_days = 8*365;          % length of simulation in days 
Tmax = Tmax_days*secpday;   % length of simulation in seconds
dt   = 3600;                % timestep in seconds
nt   = Tmax/dt;             % number of time steps in simulation
t    = dt*[0:nt];           % initilized time vector in seconds
% Time convention: t(1) = 0, t(i) = (i-1)*dt, t(nt+1) = nt*dt 

case_01p1.tdays = t/secpday;  % time vector in days
case_01p1.dz    = 10;         % box thickness [m]
case_01p1.kgas_param = 1;     % choice of air-sea flux parameterization
% integer between 1 and 6 where
% 1  Nightingale et al. (2000)
% 2  McGillis et al. (2001)
% 3  McGillis et al. (2004)
% 4  Ho et al. (2006)
% 5  Wanninkhof et al. (2009)
% 6  Wanninkhof (1992)
case_01p1.fkgas = 1;     % default is 1
% fkgas is a multiplier of kgas to explore sensitivity to uncertainty in
% kgas parameterizations
%
case_01p1.title = 'Base case';
case_01p1.legend = {' '};
% initial conditions
case_01p1.DIC_cf = 2050; % DIC counterfactual [mmol /m3]
case_01p1.alk_cf = 2300; % alk counterfactual [mmol equ /m3]
case_01p1.T = 20;        % [degC] valid range 0 - 35 degC
case_01p1.S = 35;        % valid range 20 - 40
% Typical average ocean surface temperature and salinity are 20degC and 35.
% Valid ranges for T and S reflect validity of the carbonate system 
% coefficients used in the csys routines.
%
% Set boundary conditions & forcing.
case_01p1.pCO2_air  = f_csys_alk_DIC(case_01p1.T,case_01p1.S,case_01p1.alk_cf,case_01p1.DIC_cf);  % [uatm]
case_01p1.U10  = 6.64;   % [m/s] windspeed at 10 m; global avg is 6.64 (Archer & Jacobsen 2005)
case_01p1.dalk = 100;    % alk perturbation at t=0

% FUNCTION CALLS

if 1
    % 0) call base case
    op = f_case_01p1(case_01p1);

    % plot base case
    if plot_output
        plot_case_01p1(op)
    end

    % save output to work on later
    save('output/op_case_01p1','op')
end

if 1
    % 1) Sensitivity to T and S
    T_v = [0 0 35 35];
    S_v = [20 40 20 40];
    N = length(T_v);
    
    for i=1:N
        sens_01p1_TS(i) = case_01p1; 
        sens_01p1_TS(i).title = 'Sensitivity to T and S';
        sens_01p1_TS(i).legend = {'(T,S) = (0,20)' '(T,S) = (0,40)' '(T,S) = (35,20)' '(T,S) = (35,40)'};
        sens_01p1_TS(i).T = T_v(i);
        sens_01p1_TS(i).S = S_v(i);
        sens_01p1_TS(i).pCO2_air  = f_csys_alk_DIC(T_v(i),S_v(i),sens_01p1_TS(i).alk_cf,sens_01p1_TS(i).DIC_cf);  % [uatm]
        %
        op_sens_01p1_TS(i) = f_case_01p1(sens_01p1_TS(i));
    end
    
    % save output to work on later
    save('output/op_sens_01p1_TS','op_sens_01p1_TS')

    % plot output
    if plot_output
        plot_case_01p1(op_sens_01p1_TS)        
    end
end %if       

if 1
    % 3) Sensitivity to kgas flux parameterization
    k_param = 1:5; % Not using Wanninkhof92 (option 6) because this 
    % parameterization is not recommended anymore.
    N = length(k_param);
    
    for i=1:N
        sens_01p1_k_param(i) = case_01p1;
        sens_01p1_k_param(i).title = 'Sensitivity to gas exchange parameterization';
        sens_01p1_k_param(i).legend = {'Nightingale et al. (2000)' 'McGillis et al. (2001)' ...
            'McGillis et al. (2004)' 'Ho et al. (2006)' 'Wanninkhof et al. (2009)'};
        sens_01p1_k_param(i).kgas_param = k_param(i);
        %
        op_sens_01p1_k_param(i) = f_case_01p1(sens_01p1_k_param(i));
    end
    
    % save output to work on later
    save('output/op_sens_01p1_k_param','op_sens_01p1_k_param')

    % plot output
    if plot_output
        plot_case_01p1(op_sens_01p1_k_param)   
    end
end %if

if 1
    % 4) Sensitivity to wind speed 
    U10_v = [0.5 1 1.5]*case_01p1.U10;  % -/+ 50% windspeed
    N = length(U10_v);

    for i=1:N
        sens_01p1_U10(i) = case_01p1; 
        sens_01p1_U10(i).title = 'Sensitivity to wind speed (Night2000)';
        sens_01p1_U10(i).legend = {'50% avg wind' 'avg wind' '150% avg wind'};
        sens_01p1_U10(i).U10 = U10_v(i);
        %
        op_sens_01p1_U10(i) = f_case_01p1(sens_01p1_U10(i));
    end

    % save output to work on later
    save('output/op_sens_01p1_U10','op_sens_01p1_U10')
    
    % plot output
    if plot_output
        plot_case_01p1(op_sens_01p1_U10)
    end
end %if

if 1
    % 5) Sensitivity to magnitude of kgas
    fkgas = [0.5 1 1.5];      % -/+ 50% kgas  
    N = length(fkgas);  
    
    for i=1:N
        sens_01p1_kgas(i) = case_01p1; 
        sens_01p1_kgas(i).title = 'Sensitivity to magnitude of kgas (Night2000)';
        sens_01p1_kgas(i).legend = {'0.5x kgas' 'kgas' '1.5x kgas'};
        sens_01p1_kgas(i).fkgas = fkgas(i);
        %
        op_sens_01p1_kgas(i) = f_case_01p1(sens_01p1_kgas(i));
    end
    
    % save output to work on later
    save('output/op_sens_01p1_kgas','op_sens_01p1_kgas')

    % plot output
    if plot_output
        plot_case_01p1(op_sens_01p1_kgas)
    end
end %if

if 1    
    % 6) Sensitivity to box thickness (shoud be linear where doubling 
    % the box thickness decreases the flux by half, but good to double-check)
    dz = [0.5 1 1.5]*case_01p1.dz;      % -/+ 50% dz  
    N = length(dz);  

    for i=1:N
        sens_01p1_dz(i) = case_01p1; 
        sens_01p1_dz(i).title = 'Sensitivity to box thickness';
        sens_01p1_dz(i).legend = {'0.5x dz' 'dz' '1.5x dz'};
        sens_01p1_dz(i).dz = dz(i);
        %
        op_sens_01p1_dz(i) = f_case_01p1(sens_01p1_dz(i));
    end
    % save output to work on later
    save('output/op_sens_01p1_dz','op_sens_01p1_dz')

    % plot output
    if plot_output
        plot_case_01p1(op_sens_01p1_dz)
    end
end %if

if 1    
    % 7) Sensitivity to dalk (size of alk perturbation)
    dalk = [0.5 1 1.5]*case_01p1.dalk;      % -/+ 50% dalk  
    N = length(dalk);  

    for i=1:N
        sens_01p1_dalk(i) = case_01p1; 
        sens_01p1_dalk(i).title = 'Sensitivity to Delta Alk';
        sens_01p1_dalk(i).legend = {'0.5x dalk' 'dalk' '1.5x dalk'};
        sens_01p1_dalk(i).dalk = dalk(i);
        %
        op_sens_01p1_dalk(i) = f_case_01p1(sens_01p1_dalk(i));
    end
    % save output to work on later
    save('output/op_sens_01p1_dalk','op_sens_01p1_dalk')

    % plot output
    if plot_output
        plot_case_01p1(op_sens_01p1_dalk)
    end
end %if


if 1
%----------------------------------------------------------------------
% In case 01.2 the box is not in equilibrium with the atmosphere before 
% the alkalinity perturbation at t=0.
%----------------------------------------------------------------------

% The input structure has the same elements as in case 01.1. The values
% of DIC are purturbed compared to case 01.1.

% The output structure has additional elements, because the counterfactual
% evolves in time.

% degree to which seawater is under- or oversaturated
dpCO2_v = [-200 -100 0 100 200];
N = length(dpCO2_v);
% build string for legend
for i=N:-1:1 % going backwards is a trick to preallocate the arrays
    legend_str{i} = [num2str(dpCO2_v(i)) ' uatm'];
end

for i=N:-1:1 % going backwards is a trick to preallocate the arrays
    case_01p2(i) = struct(case_01p1); % initialize structure with values from case 01.1
    case_01p2(i).DIC_cf = f_csys_alk_pCO2(case_01p2(i).T,case_01p2(i).S,case_01p2(i).alk_cf,case_01p2(i).pCO2_air+dpCO2_v(i));
    case_01p2(i).title = 'Case 01.2 Over-/undersaturated initially by';
    case_01p2(i).legend = legend_str;
    %
    op_case_01p2(i) = f_case_01p2(case_01p2(i));
end

% save output to work on later
save('output/op_case_01p2','op_case_01p2')

% call plotting function
if plot_output
    plot_case_01p2(op_case_01p2)
end

end % if
