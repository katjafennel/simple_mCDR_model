function output = f2_case_01p2(c)
%-----------
% Case 1.2
%-----------

% Here the ocean is under- or oversaturated wrt to the atmosphere before 
% the alkalinity addition at t=0. The counterfactual needs to be evolved 
% in time.

% time-stepping params
secpday = 24*3600;            % seconds per day 
Nmax = length(c.tdays)-1;     % maximum number of timesteps
dt = (c.tdays(2) - c.tdays(1))*secpday; % timestep in seconds
U10  = c.U10;                 % wind speed at 10 m

% other initilizations
Ko   = Ko_Weiss(c.T,c.S);     % [mmol /m3 /atm] solubility of CO2
% kgas [m/s] gas transfer coefficient
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
        Sc=Schmidt(T);                   
        if U10<=6    
            kgas = 0.31*U10*U10*((Sc/660)^-0.5);
        else
            kgas = 0.39*U10*U10*((Sc./660)^-0.5);
        end
end
% c.fkgas is a multiplier of kgas to explore sensitivity to uncertainty in
% kgas parameterizations; default is 1;

% initialize state vars (alk doesn't evolve in time here, but I'm anticipating
% future cases where it will; pCO2 not strictly a state var but useful to
% save as diagnostic)
% realistic member of the pair
DIC = zeros(1,Nmax+1); DIC(1) = c.DIC_cf;
alk = zeros(1,Nmax+1); alk(1) = c.alk_cf + c.dalk;
pCO2 = zeros(1,Nmax+1);      % [uatm]
% counterfactual
DIC_cf = zeros(1,Nmax+1); DIC_cf(1) = c.DIC_cf;
alk_cf = zeros(1,Nmax+1); alk_cf(1) = c.alk_cf;
pCO2_cf = zeros(1,Nmax+1);      % [uatm]

% initialize output (in addition to state vars)
Fair_sea = zeros(1,Nmax);               % air-sea flux [mmol CO2 /m2 /s] 
Fair_sea_cf = zeros(1,Nmax);            % air-sea flux [mmol CO2 /m2 /s] 
% Fair_sea > 0 is flux into the ocean (from air to sea) unlike in
% Wanninkhof et al. (2009)

% 1) Calculate integrated net air-sea flux assuming 
%    the system has reached equilibrium 
% --------------------------------------------------

% Check that difference in seawater pCO2 between the counterfactual 
% and the air are what they should be to make sure that the
% initilizations were done correctly.
%

pCO2_cf(1) = f_csys_alk_DIC(c.T,c.S,alk_cf(1),DIC_cf(1));

check = c.pCO2_air - pCO2_cf(1);
disp([' Initilization check. Is the pCO2 perturbation what it should be? dpCO2 = ' ...
        num2str(check)])

% Calculate DICdeficit/maxCDR and efficiency 

% DIC using the initial values for everything except for alk, 
% which uses the new alk (initial alk + OAE perturbation). This is the
% DIC value when air-sea gas exchange due to OAE has come to completion.
% 
DIC_oae = f_csys_alk_pCO2(c.T,c.S,alk(1),c.pCO2_air); % realistic case
DIC_oae_cf = f_csys_alk_pCO2(c.T,c.S,alk_cf(1),c.pCO2_air); % counterfactual
maxCDR = DIC_oae - DIC_oae_cf;
eta = maxCDR/c.dalk*100;  % [%]
disp([' maxCDR = ' num2str(maxCDR) '; eta = ' num2str(eta) ' %'])

% 2) Evolve both cases by time-stepping forward
% --------------------------------------------------

% time stepping
for i=1:Nmax
    % realistic
    pCO2(i) = f_csys_alk_DIC(c.T,c.S,alk(i),DIC(i)); % [uatm]
    dpCO2 = (c.pCO2_air - pCO2(i))*1.e-6; % delta pCO2  [atm]
    Fair_sea(i) = c.fkgas*kgas*Ko*dpCO2; % [mmol / m2 /s] 
    DIC(i+1) = DIC(i) + dt*Fair_sea(i)/c.dz;  
    alk(i+1) = alk(i);
    % counterfactual
    pCO2_cf(i) = f_csys_alk_DIC(c.T,c.S,alk_cf(i),DIC_cf(i)); % [uatm]
    dpCO2 = (c.pCO2_air - pCO2_cf(i))*1.e-6; % delta pCO2  [atm]
    Fair_sea_cf(i) = c.fkgas*kgas*Ko*dpCO2; % [mmol / m2 /s] 
    DIC_cf(i+1) = DIC_cf(i) + dt*Fair_sea_cf(i)/c.dz;  
    alk_cf(i+1) = alk_cf(i);
end

% 3) Calculate time to 50%, 90%, 99% of maxCDR
% --------------------------------------------------

pctOAE = [0.5 0.9 0.99];   % [50% 90% 99%]
tau = NaN*ones(size(pctOAE));
for i=1:length(tau)
    delta_DIC = DIC - DIC_cf;
    ind = find(delta_DIC >= pctOAE(i)*maxCDR); 
    if ~isempty(ind) 
        tau(i) = c.tdays(ind(1)); % [days] time until certain percentage of maxCDR is achieved
    end 
end
% calculate efficiency at t99% (note that this depends on the loop above;
% not good coding pratice)
efficiency = delta_DIC(ind(3))/maxCDR;

output.maxCDR      = maxCDR;     % same as CDRdeficit in this case
output.eta         = eta;        % isocapnic quotient
output.efficiency  = efficiency;
output.tau         = tau;
output.pCO2        = pCO2;
output.DIC         = DIC;
output.Fair_sea    = Fair_sea;
output.pCO2_cf     = pCO2_cf;
output.DIC_cf      = DIC_cf;
output.Fair_sea_cf = Fair_sea_cf;
output.tdays       = c.tdays;    % same as input; used for plotting
output.title       = c.title;    % same as input; used for plotting
output.legend      = c.legend;   % same as input; used for plotting
output.input       = c;          % save complete input parameter set

end
function [Sc]=Schmidt(T)
    A = 2073.1;     B = 125.62;     C = 3.6276;     D = 0.043219;
    Sc= A - (B.*T)+(C.*T.^2)-(D.*T.^3);
end