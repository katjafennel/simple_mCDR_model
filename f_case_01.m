function [maxCDR,efficiency,tau,tdays,pCO2,DIC,Fair_sea] = f_case_01(T,S,DIC_cf,alk_cf,pCO2_air,U10,dalk,kgas_param,fkgas)

% box thickness
dz = 10;                    % thickness [m]
%disp([' Box thickness is ' num2str(dz)])

% tau = [t_50 t_90 t_99] timescales to 50%, 90%, and 99% efficiency

% time-stepping coeffs
secpday   = 24*3600;        % seconds per day 
Tmax_days = 5*365;          % length of simulation in days 
Tmax = Tmax_days*secpday;   % length of simulation in seconds
dt   = 3600;                % timestep in seconds
Nmax = Tmax/dt;             % maximum number of timesteps
t    = dt*[0:Nmax];         % initilize time vector in seconds
% Time convention: t(1) = 0, t(i) = (i-1)*dt, t(Nmax+1) = Nmax*dt 
tdays = t/secpday;          % time vector in days

Ko   = Ko_Weiss(T,S);       % [mmol /m3 /atm] solubility of CO2
% kgas [m/s] gas transfer coefficient
switch kgas_param
    case 1 %  Nightingale et al. (2000)
        kgas = 0.333*U10 +0.222*U10*U10; 
    case 2 %  Ho et al. (2006)
        kgas = 0.266*U10*U10;            
    case 3 %  McGillis et al. (2001)
        kgas = 3.3+0.026*U10*U10*U10;    
    case 4 %  McGillis et al. (2004)
        kgas = 8.2+0.014*U10*U10*U10;    
    case 5 % Wanninkhof (1992)
        Sc=Schmidt(T);                   
        if U10<=6    
            kgas = 0.31*U10*U10*((Sc/660)^-0.5);
        else
            kgas = 0.39*U10*U10*((Sc./660)^-0.5);
        end
end
% fkgas is a multiplier of kgas to explore sensitivity to uncertainty in
% kgas parameterizations; default is fkgas = 1;

% initialize state vars (alk doesn't change here yet, but I'm anticipating
% future cases where it will; pCO2 not strictly a state var but useful to
% save as diagnostic)
DIC = zeros(1,Nmax+1); DIC(1) = DIC_cf;
alk = zeros(1,Nmax+1); alk(1) = alk_cf + dalk;
pCO2 = zeros(1,Nmax+1);      % [uatm]

% initialize output (in addition to state vars)
Fair_sea = zeros(1,Nmax);               % air-sea flux [mmol CO2 /m2 /s]   
% Fair_sea > 0 is flux into the ocean (from air to sea) unlike in
% Wanninkhof et al. (2009)
Cumul_Fair_sea = zeros(1,Nmax+1);       % cumulative Fair-sea [mmol CO2 /m2]

% 1) Calculate integrated net air-sea flux assuming 
%    the system has reached equilibrium 
% --------------------------------------------------

% Given T, S, alk, and pCO2, is the DIC returned by f_csys_alk_pCO2 the 
% same as DIC_cf? This is just a consistency check of the csys routines.
%

DIC_check = f_csys_alk_pCO2(T,S,alk_cf,pCO2_air);
err = abs(DIC_check - DIC_cf);
if err < 1E-05
    %disp(' Initilization consistent.' )
else
    disp([' ACHTUNG: Initilization inconsistent. Absolute difference: ' ...
        num2str(err)])
end

% Calculate DIC using the initial values for everything except for alk, 
% which uses the new alk (initial alk + OAE perturbation).
% 
DIC_oae = f_csys_alk_pCO2(T,S,alk(1),pCO2_air);
maxCDR = DIC_oae - DIC_cf;
efficiency = maxCDR/dalk*100;  % in %
disp([' maxCDR = ' num2str(maxCDR) ...
    '; efficiency = ' num2str(efficiency) ' %'])

% 2) Calculate integrated net air-sea flux by time-
%    stepping forward
% --------------------------------------------------

% time stepping
for i=1:Nmax
    pCO2(i) = f_csys_alk_DIC(T,S,alk(i),DIC(i)); % [uatm]
    dpCO2 = (pCO2_air - pCO2(i))*1.e-6; % delta pCO2  [atm]
    Fair_sea(i) = fkgas*kgas*Ko*dpCO2; % [mmol / m2 /s] 
    DIC(i+1) = DIC(i) + dt*Fair_sea(i)/dz;  
    alk(i+1) = alk(i);
    Cumul_Fair_sea(i+1) = Cumul_Fair_sea(i) + dt*Fair_sea(i)/dz;
end

% calculate time to 50%, 90%, 99% of maxCDR
tau = NaN*ones(1,3);
ind = find(Cumul_Fair_sea >= 0.5*maxCDR); if ~isempty(ind) tau(1) = tdays(ind(1)); end % [days] time to 50% maxCDR achieved
ind = find(Cumul_Fair_sea >= 0.9*maxCDR); if ~isempty(ind) tau(2) = tdays(ind(1)); end % [days] time to 90% maxCDR achieved
ind = find(Cumul_Fair_sea >= 0.99*maxCDR); if ~isempty(ind) tau(3) = tdays(ind(1)); end % [days] time to 99% maxCDR achieved

end
function [Sc]=Schmidt(T)
    A = 2073.1;     B = 125.62;     C = 3.6276;     D = 0.043219;
    Sc= A - (B.*T)+(C.*T.^2)-(D.*T.^3);
end