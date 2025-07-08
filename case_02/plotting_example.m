%
% example script for how to produce nice plots for case 2 runs
%

clear; close all

secpday = 24*3600;                 % seconds per day 
tmax = 400;                        % plot up to tmax (days)
fs = 18;                           % fontsize
mildBlue = [0      0.4470 0.7410]; % a pleasant blue
mildRed  = [0.8500 0.3250 0.0980]; % a pleasant red
drkYell  = [0.9290 0.6940 0.1250]; % dark yellow
purple   = [0.4940 0.1840 0.5560];
mildGrn  = [0.4660 0.6740 0.1880]; % a pleasant green
cyan     = [0.3010 0.7450 0.9330];

color_sequence = [mildBlue; purple; mildRed; drkYell; mildGrn; cyan];

label = {'Base case (k_v = 0)' 'k_v = 10^{-7} s^{-1}; \Delta_v DIC = 0' ...
    'k_v = 10^{-6} s^{-1}; \Delta_v DIC = 0' 'k_v = 10^{-5} s^{-1}; \Delta_v DIC = 0' ...
    'k_v = 10^{-6} s^{-1}; \Delta_v DIC = 100 uM' 'k_v = 10^{-6} s^{-1}; \Delta_v DIC = 200 uM'};
fn_root = {'02p2_kv0' '02_kv1em7' '02_kv1em6' '02_kv1em5' '02p3_2150' '02p3_2250'};
op_label = {'op_base' 'op_kv1em7' 'op_kv1em6' 'op_kv1em6' 'op_2150' 'op_2250'};


% initialize vectors  
fn = ['output/op_case_' fn_root{1}]; 
load(fn)
CDR_F = NaN*ones(length(label), length(op(2).t)-1+length(op(3).t)-1);
maxCDR = length(label);
tau = NaN*ones(length(label),3);

% loop over all cases
for i=1:length(label)

    disp(['Working on ' label{i}])
    fn = ['output/op_case_' fn_root{i}];
    load(fn)

    % concatenate CDR from res2 and res3
    t_CDR = [op(2).t(1:end-1) op(3).t(1:end-1)+op(2).t(end)];
    CDR_F2 = cumsum(op(2).F_int(1:end)-op(2).F_cf_int(1:end))*case_02.dt(2);
    CDR_F3 = cumsum(op(3).F_int(1:end)-op(3).F_cf_int(1:end))*case_02.dt(3) + CDR_F2(end);
    CDR_F(i,:) = [CDR_F2' CDR_F3'];
    
    % calculate time to 50%, 90%, 99% of maxCDR
    maxCDR(i) = CDR_F(i,end);
    ind = find(CDR_F(i,:) >= 0.5*maxCDR(i)); if ~isempty(ind) tau(i,1) = t_CDR(ind(1))/secpday; end % [days] time to 50% maxCDR achieved
    ind = find(CDR_F(i,:) >= 0.9*maxCDR(i)); if ~isempty(ind) tau(i,2) = t_CDR(ind(1))/secpday; end % [days] time to 90% maxCDR achieved
    ind = find(CDR_F(i,:) >= 0.99*maxCDR(i)); if ~isempty(ind) tau(i,3) = t_CDR(ind(1))/secpday; end % [days] time to 99% maxCDR achieved
    
end

figure
i=1; plot(t_CDR/secpday,CDR_F(i,:)/100,'LineWidth',2,'Color',color_sequence(i,:))
hold on
for i=2:6
    plot(t_CDR/secpday,CDR_F(i,:)/100,'LineWidth',2,'Color',color_sequence(i,:))
end
ylabel('CDR(t)/V (mmol/m3)')
xlabel('time (days)')
set(gca,'FontSize',fs)
%
Va = axis;
Va(2) = tmax;
Va(4) = 85.00;
axis(Va)

legend(label{1:6},'Location','best')
