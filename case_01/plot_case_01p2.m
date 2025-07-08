function plot_case_01p2(op)
%
% Plot output from case 01.2. Here the DIC in both realistic and counterfactual
% evolve in time as well because the seawater is over- or understaurated at the
% beginning of the simulation.
%

title_str  = op(1).title;
legend_str = op(1).legend;
td         = op(1).tdays;
Nmax       = length(op(1).tdays)-1 ;
N          = length(op);

for i=1:N
    
figure
plot(td(1:end-1),op(i).pCO2(1:end-1),td(1:end-1),op(i).pCO2_cf(1:end-1))
ylabel('pCO_2')
xlabel('time in days')
legend('real','cf')
title([title_str legend_str(i)])

figure
plot(td,op(i).DIC,td,op(i).DIC_cf)
ylabel('DIC')
xlabel('time in days')
legend('real','cf')
title([title_str legend_str(i)])

figure
plot(td(1:end-1),cumsum(op(i).Fair_sea),td(1:end-1),cumsum(op(i).Fair_sea_cf))
ylabel('Cumulative air-sea flux')
xlabel('time in days')
legend('real','cf')
title([title_str legend_str(i)])

figure
plot(td(1:end-1),op(i).Fair_sea,td(1:end-1),op(i).Fair_sea_cf)
ylabel('Air-sea flux')
xlabel('time in days')
legend('real','cf')
title([title_str legend_str(i)])

end

figure % pCO2
for i=1:N
    plot(td(1:end-1),op(i).pCO2_cf(1:end-1)-op(i).pCO2(1:end-1))
    hold on
end
ylabel('pCO_2 cf - pCO_2')
xlabel('time in days')
legend(legend_str)
title([title_str])

figure % DIC
for i=1:N
    plot(td,op(i).DIC-op(i).DIC_cf)
    hold on
end
ylabel('DIC - DIC cf')
xlabel('time in days')
legend(legend_str)
title([title_str])

    
end