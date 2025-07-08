function op = f_case_02p1(ip)

big_output = 0;

K = ip.K;
DIC_cf = ip.DIC_cf;    % DIC counterfactual [mmol /m3]
alk_cf = ip.alk_cf;    % alk counterfactual [mmol equ /m3]
dalk   = ip.dalk;      % alk perturbation at t=0 [mmol equ /m3]

res = 1;

disp([' Resolution ' num2str(res) ', dx = ' num2str(ip.dx(res)) ' m'])

% initilize time and space domains
dt = ip.dt(res);
nt = ip.nt(res);
t  = dt*[0:nt];         % initilized time vector in seconds
nx = ip.nx(res);
dx = ip.dx(res);
mp = ip.mp(res);        % midpoint of the domain
x = dx*([1:nx]-mp);     % space domain in m
dz = ip.dz;

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

disp(['Take time of time stepping loop no: ' num2str(res)])
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

% saving results in op structure
op(res).DIC_end = DIC_pi(end,:)/pi; % only save final distribution
op(res).alk_end = alk_pi(end,:)/pi; % only save final distribution
op(res).t = t;
op(res).x = x;     


% Loop over two levels of lower resolution (increase dx)
for res=2:3

    disp([' Resolution ' num2str(res) ', dx = ' num2str(ip.dx(res)) ' m'])

    % initilize time and space domains
    dt = ip.dt(res);         % timestep in seconds 
    nt = ip.nt(res);         % number of time steps in simulation
    t  = dt*[0:nt];          % initilized time vector in seconds
    nx = ip.nx(res);         % number of grid points in domain
    dx = ip.dx(res);         % delta x (m)
    mp = ip.mp(res);         % mid-point of the domain
    x = dx*([1:nx]-mp);      % domain [m]

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
    
    % resample the final distribution from the next higher res simulation (hr)
    x_hr = op(res-1).x;
    DIC_hr = op(res-1).DIC_end;
    alk_hr = op(res-1).alk_end;
    for i=1:length(x_hr)
        ind = find(x == x_hr(i));
        DIC_pi(1,ind) = DIC_hr(i)*pi;
        alk_pi(1,ind) = alk_hr(i)*pi;
    end

    disp('Take time of 2nd time stepping loop.')
    tic
    % the time stepping loop
    for i=2:nt
        
        % gas-exchange
        pCO2 = pCO2_fun(DIC_pi(i-1,:)/pi,alk_pi(i-1,:)/pi,ip.pCO2_dict);
        dpCO2 = (ip.pCO2_air - pCO2)*1.e-6; % delta pCO2  [atm]        
        Fair_sea(i,:) = ip.kgas*ip.Ko*dpCO2; % [mmol / m2 /s] 
        DIC_pi(i,:) = DIC_pi(i-1,:) + dt*pi*Fair_sea(i,:)/dz;  
    
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
        op(res).DIC = DIC_pi/pi; 
        op(res).alk = alk_pi/pi;
    end
    op(res).DIC_end = DIC_pi(end,:)/pi; % only save final distribution
    op(res).alk_end = alk_pi(end,:)/pi; % only save final distribution
    op(res).CDR = CDR;   % total mol
    op(res).t = t;
    op(res).x = x;
end

end