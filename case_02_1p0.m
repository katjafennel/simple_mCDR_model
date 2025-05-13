%
% Case 2 with dispersion
%

clear; close all

save_output = 1;

% Initilizations and preparation

secpday   = 24*3600;% seconds per day 

case_02.DIC_cf = 2050;      % DIC counterfactual [mmol /m3]
case_02.alk_cf = 2300;      % alk counterfactual [mmol equ /m3]
case_02.dalk = 100;         % desired alk perturbation at t=0 [mmol equ /3]

case_02.T = 20;             % [degC] valid range 0 - 35 degC
case_02.S = 35;             % valid range 20 - 40
case_02.pCO2_air  = f_csys_alk_DIC(case_02.T,case_02.S,case_02.alk_cf,case_02.DIC_cf);  % [uatm]
case_02.U10  = 6.64;        % [m/s] windspeed at 10 m; global avg is 6.64 (Archer & Jacobsen 2005)
case_02.Ko   = Ko_Weiss(case_02.T,case_02.S);     % [mmol /m3 /atm] solubility of CO2
case_02.kgas = 0.333*case_02.U10 +0.222*case_02.U10*case_02.U10; %  Nightingale et al. (2000)
case_02.dz = 10;            % assumed surface layer thickness
case_02.K = 80;             % horizontal dispersion coefficient in m2/s

% setting up space and time domains for three levels
case_02.dx = [10 100 1000];         % horizontal resolution in m 
case_02.nx = [600 1000 2000];       % number of horizontal grid points 
case_02.mp = case_02.nx/2;          % index of midpoint of the domain
case_02.dt = [10*10*0.5/case_02.K*0.7  50  5000];  % timestep in seconds
%            (dx*dx*0.5)/K is stability criterion threshold              
case_02.nt = [1000 4*secpday/50 2500*secpday/5000];     % number of time steps in simulation

% create pCO2 lookup table for fast calculation of pCO2
DIC_min = case_02.DIC_cf;
DIC_max = case_02.DIC_cf + case_02.dalk;
alk_min = case_02.alk_cf;
alk_max = case_02.alk_cf + case_02.dalk;
del = 1; % accuracy of the dictionary [uM]
case_02.pCO2_dict = pCO2_dictionary(case_02.T,case_02.S,DIC_min,DIC_max,alk_min,alk_max,del);


if 1  % Base case
    disp(' ')
    disp('%%% Base case')
    disp(' ')
    % Call function for timestepping
    op_case_02 = f_case_02(case_02);
    
    
    % Saving output
    if save_output
        disp('Take time of saving output.')
        tic
        % save input & output to work on later
        save('op_case_02','case_02','op_case_02') 
        toc
    end

end

if 1  % Sensitivity to dz
    disp(' ')
    disp('%%% Sensitivity to dz')
    disp(' ')
    case_02_halfdz = case_02;
    case_02_halfdz.dz = case_02.dz/2;
    case_02_halfdz.dalk = case_02.dalk*2; % so that total alk added remains unchanged (it is simply more diluted)

    % Call function for timestepping
    op_case_02_halfdz = f_case_02(case_02_halfdz);
    
    
    % Saving output
    if save_output
        disp('Take time of saving output.')
        tic
        % save input & output to work on later
        save('op_case_02_halfdz','case_02_halfdz','op_case_02_halfdz') 
        toc
    end

end

if 1  % Sensitivity to K_h (half)
    disp(' ')
    disp('%%% Sensitivity to K_h (half)')
    disp(' ')
    case_02_halfK = case_02;
    case_02_halfK.K = case_02.K/2;

    % Call function for timestepping
    op_case_02_halfK = f_case_02(case_02_halfK);
        
    % Saving output
    if save_output
        disp('Take time of saving output.')
        tic
        % save input & output to work on later
        save('op_case_02_halfK','case_02_halfK','op_case_02_halfK') 
        toc
    end
end

if 1  % Sensitivity to K_h (double)
    disp(' ')
    disp('%%% Sensitivity to K_h (double)')
    disp(' ')    
    case_02_dblK = case_02;

    K = case_02.K*2; % doubling K
    case_02_dblK.K = K;
    case_02_dblK.dt = case_02.dx.*case_02.dx*0.5/K * 0.8;  % timestep in seconds
    %         (dx*dx*0.5)/K is stability criterion threshold; 80% should be safe 
    case_02_dblK.nt = case_02.nt*2;     % number of time steps in simulation
    case_02_dblK.nx = case_02.nx*2;       % number of horizontal grid points 

    % Call function for timestepping
    op_case_02_dblK = f_case_02(case_02_dblK);
    
    % Saving output
    if save_output
        disp('Take time of saving output.')
        tic
        % save input & output to work on later
        save('op_case_02_dblK','case_02_dblK','op_case_02_dblK') 
        toc
    end

end

if 1  % Sensitivity to dalk
    disp(' ')
    disp('%%% Sensitivity to dalk')
    disp(' ')
    case_02_halfdalk = case_02;
    case_02_halfdalk.dalk = case_02.dalk/2;

    % Call function for timestepping
    op_case_02_halfdalk = f_case_02(case_02_halfdalk);
    
    
    % Saving output
    if save_output
        disp('Take time of saving output.')
        tic
        % save input & output to work on later
        save('op_case_02_halfdalk','case_02_halfdalk','op_case_02_halfdalk') 
        toc
    end

    case_02_dbldalk = case_02;
    case_02_dbldalk.dalk = case_02.dalk*2;

    % Call function for timestepping
    op_case_02_dbldalk = f_case_02(case_02_dbldalk);
    
    
    % Saving output
    if save_output
        disp('Take time of saving output.')
        tic
        % save input & output to work on later
        save('op_case_02_dbldalk','case_02_dbldalk','op_case_02_dbldalk') 
        toc
    end

end