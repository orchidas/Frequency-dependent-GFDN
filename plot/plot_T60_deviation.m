function [figH] = plot_T60_deviation(fs, couplingMatrix, tau, b, a, col, figH)
%%
% Plots the deviation in desired T60 from apprioximate paraunitarity of
% coupling matrix.
% Inputs:
% fs - sampling frequency
% couplingMatrix - nGrp x nGrp degree
% tau - delay line length in samples
% b, a - absorption filter coefficients (nGrp x filter order)
% col - colours to be used to plot
% figH - figure handle
% Outputs : 
% figH - figure handle to be returned
%%

nGrp = size(couplingMatrix, 1);
degree = size(couplingMatrix, 3);

S = svdPerBin(fft(couplingMatrix,2*degree,3));
min_dev = min(min(abs(S))).^(1/min(tau));
max_dev = max(max(abs(S))).^(1/max(tau));

Nfreq = 1024;
w = linspace(0,pi,Nfreq).';
H = zeros(Nfreq, nGrp);
T60_lims_min = zeros(Nfreq, nGrp);
T60_lims_max = zeros(Nfreq, nGrp);

for k = 1:nGrp
    H(:,k) = freqz(b(k,:), a(k,:), w);
    T60_lims_max(:,k) = get_T60_deviation(max_dev, fs, H(:,k));
    T60_lims_min(:,k) = get_T60_deviation(min_dev, fs, H(:,k));

    loglog((w/pi) * (fs/2),abs(H(:,k) - T60_lims_min(:,k)), 'LineWidth', ...
        1.2, 'Color', col(k,:));hold on;grid on;
    loglog((w/pi) * (fs/2),abs(T60_lims_max(:,k) - H(:,k)), 'LineWidth', ...
        1.2, 'Color', col(k,:));hold on;grid on;
end

hold off;
ylabel('T60 error (s)');
xlabel('Frequency (Hz)');
xlim([20 20000]);
set(gca, 'FontUnits','points', 'FontWeight','normal', 'FontSize',8, 'FontName','Times');


%% helper functions

function T60_new = get_T60_deviation(sigma, fs, T60_des)
    c = log(0.001)/fs;
    T60_new = c./(log(sigma) + c./T60_des);
end

end

