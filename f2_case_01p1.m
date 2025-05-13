function output = f2_case_01p1(c)
%-----------
% Case 1.1
%-----------

% Here the ocean and atmosphere are in equilibrium before the 
% alkalinity addition at t=0 and the counterfactual is in equilibrium 
% at all times (no need to step this forward in time).

% time-stepping params
secpday = 24*3600;            % seconds per day 
Nmax = length(c.tdays)-1;     % maximum number of timesteps
dt = (c.tdays(2) - c.tdays(1))*secpday; % timestep in seconds

% other initilizations
Ko   = Ko_Weiss(c.T,c.S);     % [mmol /m3 /atm] solubility of CO2
% kgas [m/s] gas transfer coefficient
U10  = c.U10;
switch c.kgas_param
    case 1 %  Nightingale et al. (2000)
        kgas = 0.333*U10 +0.222*U10*U10; 
    case 2 %  McGillis et al. (2001)
        kgas = 3.3+0.026*U10*U10*U10; 
    case 3 %  McGillis et al. (2004)
        kgas = 8.2+0.014*U10*U10*U10; 
    case 4 %  Ho et al. (2006)
        kgas = 0.266*U10*U10;            
    case 5 % Wanninkhof et al. (2009)
        kgas = 0.24*U10*U10; 
    case 6 % Wanninkhof (1992)
        Sc=Schmidt(c.T);                   
        if U10<=6    
            kgas = 0.31*U10*U10*((Sc/660)^-0.5);
        else
            kgas = 0.39*U10*U10*((Sc./660)^-0.5);
        end
end
% c.fkgas is a multiplier of kgas to explore sensitivity to uncertainty in
% kgas parameterizations; default is 1;

% initialize state vars (alk doesn't change here yet, but I'm anticipating
% future cases where it will; pCO2 not strictly a state var but useful to
% save as diagnostic)
DIC = zeros(1,Nmax+1); DIC(1) = c.DIC_cf;
alk = zeros(1,Nmax+1); alk(1) = c.alk_cf + c.dalk;
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

DIC_check = f_csys_alk_pCO2(c.T,c.S,c.alk_cf,c.pCO2_air);
err = abs(DIC_check - c.DIC_cf);
if err < 1E-05
    %disp(' Initilization consistent.' )
else
    disp([' ACHTUNG: Initilization inconsistent. Absolute difference: ' ...
        num2str(err)])
end

% Calculate DIC using the initial values for everything except for alk, 
% which uses the new alk (initial alk + OAE perturbation). This is the
% DIC value when air-sea gas exchange due to OAE has come to completion.
% 
DIC_oae = f_csys_alk_pCO2(c.T,c.S,alk(1),c.pCO2_air);
maxCDR = DIC_oae - c.DIC_cf;
eta = maxCDR/c.dalk*100;  % [%]
disp([' maxCDR = ' num2str(maxCDR) '; eta = ' num2str(eta) ' %'])

% 2) Calculate integrated net air-sea flux by time-
%    stepping forward
% --------------------------------------------------

% time stepping
for i=1:Nmax
    pCO2(i) = f_csys_alk_DIC(c.T,c.S,alk(i),DIC(i)); % [uatm]
    dpCO2 = (c.pCO2_air - pCO2(i))*1.e-6; % delta pCO2  [atm]
    Fair_sea(i) = c.fkgas*kgas*Ko*dpCO2; % [mmol / m2 /s] 
    DIC(i+1) = DIC(i) + dt*Fair_sea(i)/c.dz;  
    alk(i+1) = alk(i);
    Cumul_Fair_sea(i+1) = Cumul_Fair_sea(i) + dt*Fair_sea(i)/c.dz;
end

% calculate time to 50%, 90%, 99% of maxCDR
tau = NaN*ones(1,3);
ind = find(Cumul_Fair_sea >= 0.5*maxCDR); if ~isempty(ind) tau(1) = c.tdays(ind(1)); end % [days] time to 50% maxCDR achieved
ind = find(Cumul_Fair_sea >= 0.9*maxCDR); if ~isempty(ind) tau(2) = c.tdays(ind(1)); end % [days] time to 90% maxCDR achieved
ind = find(Cumul_Fair_sea >= 0.99*maxCDR); if ~isempty(ind) tau(3) = c.tdays(ind(1)); end % [days] time to 99% maxCDR achieved

% calculate efficiency at t99%
efficiency = Cumul_Fair_sea(ind(3))/maxCDR;

output.maxCDR      = maxCDR;
output.eta         = eta;
output.efficiency  = efficiency; % efficiency at t99%
output.tau         = tau;
output.pCO2        = pCO2;
output.DIC         = DIC;
output.Fair_sea    = Fair_sea;
output.tdays       = c.tdays;    % same as input; used for plotting
output.title       = c.title;    % same as input; used for plotting
output.legend      = c.legend;   % same as input; used for plotting
output.input       = c;          % save complete input parameter set

end
function [Sc]=Schmidt(T)
    A = 2073.1;     B = 125.62;     C = 3.6276;     D = 0.043219;
    Sc= A - (B.*T)+(C.*T.^2)-(D.*T.^3);
end