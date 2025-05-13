%
% Case 2 with dispersion and deep ocean reservoir
%
% Output naming convention:
% case_02.2 assumes no vertical gradient in DIC (i.e. DIC_deep = DIC_cf)
% case_02.3 has vertical DIC gradient (i.e. DIC_deep > DIC_cf_surface)
% case_02.4 has over- or undersaturated surface DIC

clear; close all

save_output = 1;

% Initilizations and preparation

secpday   = 24*3600;% seconds per day 

case_02.DIC_cf = 2050;      % DIC counterfactual [mmol /m3]
case_02.DIC_deep = 2050;    % deep DIC (same for real and cf) [mmol /m3]
case_02.alk_cf = 2300;      % alk counterfactual [mmol equ /m3]
case_02.dalk = 100;         % desired alk perturbation at t=0 [mmol equ /3]

case_02.T = 20;             % [degC] valid range 0 - 35 degC
case_02.S = 35;             % valid range 20 - 40
case_02.pCO2_air  = f_csys_alk_DIC(case_02.T,case_02.S,case_02.alk_cf,case_02.DIC_cf);  % [uatm]
case_02.U10  = 6.64;        % [m/s] windspeed at 10 m; global avg is 6.64 (Archer & Jacobsen 2005)
case_02.Ko   = Ko_Weiss(case_02.T,case_02.S);     % [mmol /m3 /atm] solubility of CO2
case_02.kgas = 0.333*case_02.U10 +0.222*case_02.U10*case_02.U10; %  Nightingale et al. (2000)
case_02.dz = 10;            % assumed surface layer thickness
case_02.Kh = 80;            % horizontal dispersion coefficient in m2/s
case_02.Kv = 0;             % vertical dispersion coefficient in m2/s

% setting up space and time domains for three levels
case_02.dx = [10 100 1000];         % horizontal resolution in m 
case_02.nx = [600 1000 2000];       % number of horizontal grid points 
case_02.mp = case_02.nx/2;          % index of midpoint of the domain
case_02.dt = [10*10*0.5/case_02.Kh*0.7  50  5000];  % timestep in seconds
%            (dx*dx*0.5)/K is stability criterion threshold              
case_02.nt = [1000 4*secpday/50 1000*secpday/5000];     % number of time steps in simulation

% create pCO2 lookup table for fast calculation of pCO2
DIC_min = case_02.DIC_cf;
DIC_max = max(case_02.DIC_cf + case_02.dalk,case_02.DIC_deep);
alk_min = case_02.alk_cf;
alk_max = case_02.alk_cf + case_02.dalk;
del = 1; % accuracy of the dictionary [uM]
case_02.pCO2_dict = pCO2_dictionary(case_02.T,case_02.S,DIC_min,DIC_max,alk_min,alk_max,del);


if 1  % Base case (kv = 0; no vertical DIC gradient; in equilibrium initially)
    disp(' ')
    disp('%%% Base case')
    disp(' ')

    % Call function for timestepping
    op_case_02p2_kv0 = f_case_02p2(case_02);
    
    
    % Saving output
    if save_output
        disp(' ')
        disp('Take time of saving output.')
        tic
        % save input & output to work on later
        %save('op_case_02p2_base','case_02','op') 
        save('op_case_02p2_kv0','case_02','op_case_02p2_kv0') 
        toc
        disp(' ')
    end

end

if 0  % sensitivity to K_v
    disp(' ')
    disp('%%% Sensitivity to Kv')
    disp(' ')


    case_02_kv1em7 = case_02;
    case_02_kv1em7.Kv = 0.0000001;

    % Call function for timestepping
    op_case_02_kv1em7 = f_case_02p2(case_02_kv1em7);

    % Saving output
    if save_output
        disp(' ')
        disp('Take time of saving output.')
        tic
        % save input & output to work on later
        save('op_case_02_kv1em7','case_02_kv1em7','op_case_02_kv1em7') 
        toc
    end


    % case_02_kv2em7 = case_02;
    % case_02_kv2em7.Kv = 0.0000002;
    % 
    % % Call function for timestepping
    % op = f_case_02p2(case_02_kv2em7);
    % 
    % % Saving output
    % if save_output
    %     disp(' ')
    %     disp('Take time of saving output.')
    %     tic
    %     % save input & output to work on later
    %     save('op_case_02p2_kv2em7','case_02_kv2em7','op') 
    %     toc
    % end


    case_02_kv1em6 = case_02;
    case_02_kv1em6.Kv = 0.000001;

    % Call function for timestepping
    op_case_02_kv1em6 = f_case_02p2(case_02_kv1em6);

    % Saving output
    if save_output
        disp(' ')
        disp('Take time of saving output.')
        tic
        % save input & output to work on later
        save('op_case_02_kv1em6','case_02_kv1em6','op_case_02_kv1em6') 
        toc
    end

    case_02_kv1em5 = case_02;
    case_02_kv1em5.Kv = 0.00001;

    % Call function for timestepping
    op_case_02_kv1em5 = f_case_02p2(case_02_kv1em5);


    % Saving output
    if save_output
        disp(' ')
        disp('Take time of saving output.')
        tic
        % save input & output to work on later
        save('op_case_02_kv1em5','case_02_kv1em5','op_case_02_kv1em5') 
        toc
        disp(' ')
    end

end

if 0  % sensitivity to DIC_deep
    disp(' ')
    disp('%%% Sensitivity to DIC_deep')
    disp(' ')

    case_02p3_2150 = case_02;
    case_02p3_2150.Kv = 0.000001; 
    case_02p3_2150.DIC_deep = 2150;

    % Call function for timestepping
    op_case_02p3_2150 = f_case_02p2(case_02p3_2150);

    % Saving output
    if save_output
        disp(' ')
        disp('Take time of saving output.')
        tic
        % save input & output to work on later
        save('op_case_02p3_2150','case_02p3_2150','op_case_02p3_2150') 
        toc
    end


    case_02p3_2250 = case_02;
    case_02p3_2250.Kv = 0.000001; 
    case_02p3_2250.DIC_deep = 2250;

    % create pCO2 lookup table for fast calculation of pCO2
    DIC_min = case_02.DIC_cf;
    DIC_max = max(case_02.DIC_cf + case_02.dalk,case_02.DIC_deep);
    alk_min = case_02.alk_cf;
    alk_max = case_02.alk_cf + case_02.dalk;
    del = 1; % accuracy of the dictionary [uM]
    case_02p3_2250.pCO2_dict = pCO2_dictionary(case_02.T,case_02.S,DIC_min,DIC_max,alk_min,alk_max,del);

    % Call function for timestepping
    op_case_02p3_2250 = f_case_02p2(case_02p3_2250);

    % Saving output
    if save_output
        disp(' ')
        disp('Take time of saving output.')
        tic
        % save input & output to work on later
        save('op_case_02p3_2250','case_02p3_2250','op_case_02p3_2250') 
        toc
    end

end


if 0  % sensitivity to over- and undersaturation at the surface
    disp(' ')
    disp('%%% Sensitivity to over- and undersaturation at the surface')
    disp(' ')
 
    
    case_02p4 = case_02;
    case_02p4.Kv = 0.000001; 
    case_02p4.DIC_cf = 2150; % over by 100 uM
    case_02p4.DIC_deep = 2250;

    % create pCO2 lookup table for fast calculation of pCO2
    DIC_min = case_02.DIC_cf;
    DIC_max = max(case_02.DIC_cf + case_02.dalk,case_02.DIC_deep);
    alk_min = case_02.alk_cf;
    alk_max = case_02.alk_cf + case_02.dalk;
    del = 1; % accuracy of the dictionary [uM]
    case_02p4.pCO2_dict = pCO2_dictionary(case_02.T,case_02.S,DIC_min,DIC_max,alk_min,alk_max,del);

    % Call function for timestepping
    op = f_case_02p2(case_02p4);

    % Saving output
    if save_output
        disp(' ')
        disp('Take time of saving output.')
        tic
        % save input & output to work on later
        save('op_case_02p4_2250_2150','case_02p4','op') 
        toc
    end


    case_02p4 = case_02;
    case_02p4.Kv = 0.000001; 
    case_02p4.DIC_cf = 1950; % under by 100 uM
    case_02p4.DIC_deep = 2050;

    % create pCO2 lookup table for fast calculation of pCO2
    DIC_min = case_02p4.DIC_cf;
    DIC_max = max(case_02p4.DIC_cf + case_02p4.dalk,case_02p4.DIC_deep);
    alk_min = case_02.alk_cf;
    alk_max = case_02.alk_cf + case_02.dalk;
    del = 1; % accuracy of the dictionary [uM]
    case_02p4.pCO2_dict = pCO2_dictionary(case_02.T,case_02.S,DIC_min,DIC_max,alk_min,alk_max,del);

    % Call function for timestepping
    op = f_case_02p2(case_02p4);

    % Saving output
    if save_output
        disp(' ')
        disp('Take time of saving output.')
        tic
        % save input & output to work on later
        save('op_case_02p4_2050_1950','case_02p4','op') 
        toc
    end



end