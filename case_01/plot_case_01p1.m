function plot_case_01p1(op)
%
% Plotting output for case 01.1. Since the counterfactual is timeinvariant
% there is no need to plot it. Just plotting the realistic case (with the 
% oae perturbation).

title_str  = op(1).title;
legend_str = op(1).legend;
td         = op(1).tdays;
Nmax       = length(op(1).tdays)-1 ;
N          = length(op);

if N>1
    
    figure
    pCO2 = reshape(cell2mat({op.pCO2}),Nmax+1,N);
    plot(td(1:end-1),pCO2(1:end-1,:))
    ylabel('pCO_2')
    xlabel('time in days')
    title(title_str)
    legend(legend_str,'Location','best')

    figure
    DIC = reshape(cell2mat({op.DIC}),Nmax+1,N);
    plot(td,DIC)
    ylabel('DIC')
    xlabel('time in days')
    title(title_str)
    legend(legend_str,'Location','best')

    figure
    F = reshape(cell2mat({op.Fair_sea}),Nmax,N);
    plot(td(1:end-1),cumsum(F))
    ylabel('Cumulative air-sea flux')
    xlabel('time in days')
    title(title_str)
    legend(legend_str,'Location','best')

    figure
    plot(td(1:end-1),F)
    ylabel('Air-sea flux')
    xlabel('time in days')
    title(title_str)
    legend(legend_str,'Location','best')
    
else

    figure
    plot(td(1:end-1),op.pCO2(:,1:end-1))
    ylabel('pCO_2')
    xlabel('time in days')
    title(title_str)

    figure
    plot(td,op.DIC)
    ylabel('DIC')
    xlabel('time in days')
    title(title_str)

    figure
    plot(td(1:end-1),cumsum(op.Fair_sea))
    ylabel('Cumulative air-sea flux')
    xlabel('time in days')
    title(title_str)

    figure
    plot(td(1:end-1),op.Fair_sea)
    ylabel('Air-sea flux')
    xlabel('time in days')
    title(title_str)
    
end