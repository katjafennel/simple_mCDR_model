%
% Case 2 with dispersion
%

% Initilizations and preparation

clear; close all

big_output = 0;     % yes=1/no=0
save_output = 1;    % yes=1/no=0

secpday   = 24*3600;% seconds per day 

DIC_cf = 2050;      % DIC counterfactual [mmol /m3]
alk_cf = 2300;      % alk counterfactual [mmol equ /m3]
dalk = 100;         % desired alk perturbation at t=0 [mmol equ /3]

T = 20;             % [degC] valid range 0 - 35 degC
S = 35;             % valid range 20 - 40
pCO2_air  = f_csys_alk_DIC(T,S,alk_cf,DIC_cf);  % [uatm]
U10  = 6.64;        % [m/s] windspeed at 10 m; global avg is 6.64 (Archer & Jacobsen 2005)
Ko   = Ko_Weiss(T,S);     % [mmol /m3 /atm] solubility of CO2
kgas = 0.333*U10 +0.222*U10*U10; %  Nightingale et al. (2000)
dz = 10;            % assumed surface layer thickness
K = 80;             % dispersion coefficient in m2/s

% level 1 (highest resolution)
% setting up space and time domain
dx = 10;            % horizontal resolution in m
nx = 600;           % number of horizontal grid points 
mp = nx/2;          % index of midpoint of the domain
x = dx*([1:nx]-mp);	% space domain in m
dt = dx*dx*0.5/K *0.7;  % timestep in seconds
% dt = 70% of (dx*dx*0.5)/K (70% of stability criterion threshold)                
nt = 1000;          % number of time steps in simulation
t = dt*[0:nt];      % initilized time vector in seconds

% create pCO2 lookup table
DIC_min = DIC_cf;
DIC_max = DIC_cf + dalk;
alk_min = alk_cf;
alk_max = alk_cf + dalk;
del = 1; % accuracy of the dictionary [uM]
pCO2_dict = pCO2_dictionary(T,S,DIC_min,DIC_max,alk_min,alk_max,del);

% Start timestepping

disp([' Resolution 1 (highest), dx = ' num2str(dx) ' m'])

d = K*dt/dx/dx;    % non-dimensional dispersion number 
% test stability (d has to be less than 0.5)
if d>= 0.5
    disp( 'ACHTUNG: Diffusivity stability criterion is violated. ')
    disp(['         dispersion number = ' num2str(d) ' but must be <0.5'])
    disp( '         Decrease diffusivity or time step. Or increase dx.')
    error('ERROR')
end

% initialize tracer arrays (polar coordinate; integral over 180 degrees)
DIC_pi = pi*DIC_cf*ones(nt,nx);  % 1st index is time, 2nd index is space (distance from mp)   
alk_pi = pi*alk_cf*ones(nt,nx);
% add initial alkalinity perturbation
alk_pi(1,mp) = alk_pi(1,mp) + pi*dalk; 

disp('Take time of 1st time stepping loop.')
tic
% the time stepping loop
for i=2:nt
    
    % first calculate upstream and downstream diff terms for both tracers
    DIC_diff_us = d*(-DIC_pi(i-1,2:nx-1)+DIC_pi(i-1,1:nx-2));
    DIC_diff_ds = d*(DIC_pi(i-1,3:nx)-DIC_pi(i-1,2:nx-1));
    alk_diff_us = d*(-alk_pi(i-1,2:nx-1)+alk_pi(i-1,1:nx-2));
    alk_diff_ds = d*(alk_pi(i-1,3:nx)-alk_pi(i-1,2:nx-1));
    % apply dispersion to tracer arrays
    DIC_pi(i,2:nx-1) = DIC_pi(i-1,2:nx-1) + DIC_diff_us + DIC_diff_ds;
    alk_pi(i,2:nx-1) = alk_pi(i-1,2:nx-1) + alk_diff_us + alk_diff_ds;
%   this is equal to:    
%   C(2:nx-1) = Cold(2:nx-1) + d*(-Cold(2:nx-1)+Cold(1:nx-2)+Cold(3:nx)-Cold(2:nx-1));
     
end		% end of time stepping loop
toc 

% saving results in structure; index 1 indicates highest resolution level
op_case_02(1).DIC_end = DIC_pi(end,:)/pi; % only save final distribution
op_case_02(1).alk_end = alk_pi(end,:)/pi; % only save final distribution
op_case_02(1).t = t;
op_case_02(1).dt = dt;
op_case_02(1).nt = nt;
op_case_02(1).x = x;
op_case_02(1).dx = dx;
op_case_02(1).dz = dz;
op_case_02(1).mp = mp;
op_case_02(1).K = K;
op_case_02(1).DIC_cf = DIC_cf;      % DIC counterfactual [mmol /m3]
op_case_02(1).alk_cf = alk_cf;      % alk counterfactual [mmol equ /m3]
op_case_02(1).dalk = dalk;          % alk perturbation at t=0 [mmol equ /m3]


% Lower resolution (increase dx)

dx = 100;           % horizontal resolution in m
nx = 1000;          % number of horizontal grid points 
mp = nx/2;          % index of midpoint of the domain
x = dx*([1:nx]-mp);	% space domain in m
Tmax_days = 4;      % length of simulation in days 
Tmax = Tmax_days*secpday;   % length of simulation in seconds
dt = 50;            % timestep in seconds (dx*dx*0.5/K *0.8)
nt = Tmax/dt;       % number of time steps in simulation
t = dt*[0:nt];      % initilized time vector in seconds

disp([' Resolution 2 (lower), dx = ' num2str(dx) ' m'])
d = K*dt/dx/dx;    % non-dimensional dispersion number 
% test stability (d has to be less than 0.5)
if d>= 0.5
    disp( 'ACHTUNG: Diffusivity stability criterion is violated. ')
    disp(['         dispersion number = ' num2str(d) ' but must be <0.5'])
    disp( '         Decrease diffusivity or time step. Or increase dx.')
    error('ERROR')
end

% initialize tracer arrays
DIC_pi = pi*DIC_cf*ones(nt,nx);  % 1st index is time, 2nd index is space   
alk_pi = pi*alk_cf*ones(nt,nx);
% initilize diagnostics
Fair_sea = zeros(nt,nx);               % air-sea flux [mmol CO2 /m2 /s]   
% Fair_sea > 0 is flux into the ocean (from air to sea) unlike in
% Wanninkhof et al. (2009)

% resample the final distribution from the high-res (hr) simulation 
x_hr = op_case_02(1).x;
DIC_hr = op_case_02(1).DIC_end;
alk_hr = op_case_02(1).alk_end;
for i=1:length(x_hr)
    ind = find(x == x_hr(i));
    DIC_pi(1,ind) = DIC_hr(i)*pi;
    alk_pi(1,ind) = alk_hr(i)*pi;
end

disp('Take time of 2nd time stepping loop.')
tic
% the time stepping loop
for i=2:nt
    
    if 1
    % gas-exchange
    pCO2 = pCO2_fun(DIC_pi(i-1,:)/pi,alk_pi(i-1,:)/pi,pCO2_dict);
    dpCO2 = (pCO2_air - pCO2)*1.e-6; % delta pCO2  [atm]        
    Fair_sea(i,:) = kgas*Ko*dpCO2; % [mmol / m2 /s] 
    DIC_pi(i,:) = DIC_pi(i-1,:) + dt*pi*Fair_sea(i,:)/dz;  
    end

    % dispersion
    % first calculate upstream and downstream diff terms for both tracers
    DIC_diff_us = d*(-DIC_pi(i-1,2:nx-1)+DIC_pi(i-1,1:nx-2));
    DIC_diff_ds = d*(DIC_pi(i-1,3:nx)-DIC_pi(i-1,2:nx-1));
    alk_diff_us = d*(-alk_pi(i-1,2:nx-1)+alk_pi(i-1,1:nx-2));
    alk_diff_ds = d*(alk_pi(i-1,3:nx)-alk_pi(i-1,2:nx-1));
    % apply dispersion to tracer arrays
    DIC_pi(i,2:nx-1) = DIC_pi(i,2:nx-1) + DIC_diff_us + DIC_diff_ds;
    alk_pi(i,2:nx-1) = alk_pi(i-1,2:nx-1) + alk_diff_us + alk_diff_ds;
%   this is equal to:    
%   C(2:nx-1) = Cold(2:nx-1) + d*(-Cold(2:nx-1)+Cold(1:nx-2)+Cold(3:nx)-Cold(2:nx-1));
     
end		% end of time stepping loop
toc

% calculate CDR(t) in mol
CDR = sum(DIC_pi/pi-DIC_cf,2)*dx*dz;

% index 2 indicates 2nd resolution level (coarser)
if big_output
    op_case_02(2).DIC = DIC_pi/pi; 
    op_case_02(2).alk = alk_pi/pi;
end
op_case_02(2).DIC_end = DIC_pi(end,:)/pi; % only save final distribution
op_case_02(2).alk_end = alk_pi(end,:)/pi; % only save final distribution
op_case_02(2).CDR = CDR;   % total mol
op_case_02(2).t = t;
op_case_02(2).dt = dt;
op_case_02(2).nt = nt;
op_case_02(2).x = x;
op_case_02(2).dx = dx;
op_case_02(2).dz = dz;
op_case_02(2).mp = mp;
op_case_02(2).K = K;
op_case_02(2).DIC_cf = DIC_cf;      % DIC counterfactual [mmol /m3]
op_case_02(2).alk_cf = alk_cf;      % alk counterfactual [mmol equ /m3]


% Lowest resolution (increase dx one more time)

dx = 1000;          % horizontal resolution in m
nx = 2000;          % number of horizontal grid points 
mp = nx/2;          % index of midpoint of the domain
x = dx*([1:nx]-mp);	% space domain in m
Tmax_days = 1000;   % length of simulation in days 
Tmax = Tmax_days*secpday;   % length of simulation in seconds
K = 80;             % dispersion coefficient in m2/s applied in 1D (dispersion in 2D is K/pi)
dt = 5000;          % timestep in seconds (dx*dx*0.5/K *0.8)
nt = Tmax/dt;       % number of time steps in simulation
t = dt*[0:nt];      % initilized time vector in seconds

disp([' Resolution 3 (lowest), dx = ' num2str(dx) ' m'])
d = K*dt/dx/dx;    % non-dimensional dispersion number 
% test stability (d has to be less than 0.5)
if d>= 0.5
    disp( 'ACHTUNG: Diffusivity stability criterion is violated. ')
    disp(['         dispersion number = ' num2str(d) ' but must be <0.5'])
    disp( '         Decrease diffusivity or time step. Or increase dx.')
    error('ERROR')
end

% initialize tracer arrays
DIC_pi = pi*DIC_cf*ones(nt,nx);  % 1st index is time, 2nd index is space   
alk_pi = pi*alk_cf*ones(nt,nx);
% initilize diagnostics
pCO2 = zeros(nt,nx);  
Fair_sea = zeros(nt,nx);               % air-sea flux [mmol CO2 /m2 /s]   
% Fair_sea > 0 is flux into the ocean (from air to sea) unlike in
% Wanninkhof et al. (2009)

% resample the final distribution from the higher-res (hr) simulation 
x_hr = op_case_02(2).x;
DIC_hr = op_case_02(2).DIC_end; 
alk_hr = op_case_02(2).alk_end;
for i=1:length(x_hr)
    ind = find(x == x_hr(i));
    DIC_pi(1,ind) = DIC_hr(i)*pi;
    alk_pi(1,ind) = alk_hr(i)*pi;
end

disp('Take time of 3nd time stepping loop.')
tic
% the time stepping loop
for i=2:nt

    if 1
    % gas-exchange
    pCO2 = pCO2_fun(DIC_pi(i-1,:)/pi,alk_pi(i-1,:)/pi,pCO2_dict);
    dpCO2 = (pCO2_air - pCO2)*1.e-6; % delta pCO2  [atm]        
    Fair_sea(i,:) = kgas*Ko*dpCO2; % [mmol / m2 /s] 
    DIC_pi(i,:) = DIC_pi(i-1,:) + dt*pi*Fair_sea(i,:)/dz;  
    end

    % dispersion
    % first calculate upstream and downstream diff terms for both tracers
    DIC_diff_us = d*(-DIC_pi(i-1,2:nx-1)+DIC_pi(i-1,1:nx-2));
    DIC_diff_ds = d*(DIC_pi(i-1,3:nx)-DIC_pi(i-1,2:nx-1));
    alk_diff_us = d*(-alk_pi(i-1,2:nx-1)+alk_pi(i-1,1:nx-2));
    alk_diff_ds = d*(alk_pi(i-1,3:nx)-alk_pi(i-1,2:nx-1));
    % apply dispersion to tracer arrays
    DIC_pi(i,2:nx-1) = DIC_pi(i,2:nx-1) + DIC_diff_us + DIC_diff_ds;
    alk_pi(i,2:nx-1) = alk_pi(i-1,2:nx-1) + alk_diff_us + alk_diff_ds;
%   this is equal to:    
%   C(2:nx-1) = Cold(2:nx-1) + d*(-Cold(2:nx-1)+Cold(1:nx-2)+Cold(3:nx)-Cold(2:nx-1));
     
end		% end of time stepping loop
toc

% calculate CDR(t) in mol
CDR = sum(DIC_pi/pi-DIC_cf,2)*dx*dz;

% index 3 indicates 3rd resolution level (coarsest)
if big_output
    op_case_02(3).DIC = DIC_pi/pi; 
    op_case_02(3).alk = alk_pi/pi;
end
op_case_02(3).DIC_end = DIC_pi(end,:)/pi; % only save final distribution
op_case_02(3).alk_end = alk_pi(end,:)/pi; % only save final distribution
op_case_02(3).CDR = CDR;   % total mol
op_case_02(3).t = t;
op_case_02(3).dt = dt;
op_case_02(3).nt = nt;
op_case_02(3).x = x;
op_case_02(3).dx = dx;
op_case_02(3).dz = dz;
op_case_02(3).mp = mp;
op_case_02(3).K = K;
op_case_02(3).DIC_cf = DIC_cf;      % DIC counterfactual [mmol /m3]
op_case_02(3).alk_cf = alk_cf;      % alk counterfactual [mmol equ /m3]

if save_output
    % save output to work on later
    save('op_case_02','op_case_02')
end

