function [figH] = plot_desired_T60(fs, nGrp, b, a,figH)
%%
% Plot the desired T60 response
% Inputs
% b,a - T60 filter coefficients (IIR) - nGrp x K
% fs - sampling frequency
% nGrp - number of groups in the GFDN
% figH - figure handle

Nfreq = 1024;
w = linspace(0,pi,Nfreq).';
H = zeros(Nfreq, nGrp);
legendStr = {};

for k = 1:nGrp
    H(:,k) = freqz(b(k,:), a(k,:), w);
    loglog((w/pi) * (fs/2),abs(H(:,k)), 'LineWidth', 1.2);hold on;grid on;
    legendStr{k} = ['Room ', num2str(k)];
end

hold off;
legend(legendStr, 'Location', 'southwest');
ylabel('T60 desired (s)');
yticks([0 0.1 0.5 1 2 3 4]);
xlabel('Frequency (Hz)');
xlim([20 20000]); ylim([0,4]);
set(gca, 'FontUnits','points', 'FontWeight','normal', 'FontSize',8, 'FontName','Times');
end

