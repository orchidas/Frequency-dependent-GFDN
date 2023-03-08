function [] = plot_fdn_attenuation_filters(nDel, nSize, fs, b, a, col, ls)
%%
% Plot attenuation filter response in the FDN
% Inputs:
% nDel - total number of delay lines
% nSize - number of delay lines in each group
% fs - sampling frequency
% b, a - attenuation filter coefficients (nDel x filter order)
% col - colors used for plotting
% ls - linestyles used for plotting
%%

Nfreq = 1024;
w = linspace(0,pi,Nfreq).';
plt = gobjects(1,nDel);
legendStr = {};

assert (sum(nSize) == nDel, 'Number of delay lines in each group does not add up to total number of delay lines');

for i = 1:nDel
    G = freqz(b(i,:),a(i,:),w);  
    plt(i) = semilogx((w/pi) * (fs/2),20*log10(abs(G)), 'Color',col(ceil(i/nSize(1)),:), ...
    'LineWidth', 1.2, 'LineStyle',ls{ceil(i/nSize(1))});
    grid on;hold on;
    
    
    if mod(i, nSize(1)) == 0
        nGrp = i / nSize(1);
        legendStr{nGrp} = ['Room ', num2str(nGrp)];
    end
end

hold off;
ylabel('T60 Filter gain (dB)');
xlabel('Frequency (Hz)');
legend([plt(1), plt(nSize(1)+1), plt(2*nSize(1)+1)],legendStr,...
    'Location','southwest','FontSize',5);
xlim([20 30000]);ylim([-8 0]);
set(gca, 'FontUnits','points', 'FontWeight','normal', 'FontSize',8, 'FontName','Times');

