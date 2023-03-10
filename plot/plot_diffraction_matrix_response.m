function [plotAxes,plotHandles] = plot_diffraction_matrix_response(nRooms,fs,c,aperture,varargin)
%% Frequency response of designed diffraction-based coupling matrices
% INPUTS:
% nRooms - number of rooms (int), 
% fs - sampling frequency (Hz)
% c = speed of sound, (m/s)
% aperture - size (can be a vector)
% OPTIONAL:
% areas - areas of coupled rooms
% absorption_coefs - absorption coefficient of each room
% OUTPUTS:
% plotAxes, plotHandles - axes and handles of plot
    
addpath(genpath('plot/.'));
p = inputParser;
p.KeepUnmatched = true;
addParameter(p, 'areas', ones(length(aperture),1));
addParameter(p, 'absorption_coeffs', []);
addParameter(p, 'save_flag', 0);

parse(p,varargin{:});
area = p.Results.areas;
absorp_coef = p.Results.absorption_coeffs;
save_flag = p.Results.save_flag;

nbins = 2^10;
w = pi*((0:nbins-1).'/nbins);
lgdstr = {};

colors =  [[0    0.4470    0.7410];
        [0.8500    0.3250    0.0980];
        [0.9290    0.6940    0.1250];
        [0.4940    0.1840    0.5560];
        [0.4660    0.6740    0.1880];
        [0.3010    0.7450    0.9330];
        [0.6350    0.0780    0.1840]];


for k = 1:length(aperture)

     couplingMatrixTD = two_room_coupling_matrix(nRooms, aperture(k), area, absorp_coef, fs, c, 'exact', 'asymmetry', 0);
     couplingMatrixFD = fft(couplingMatrixTD, nbins*2, 3);
     
     if k < length(aperture)
        [plotAxes, plotHandles] = plotImpulseResponseMatrixStem((w/pi)*(fs/2), ...
            20*log10(abs(couplingMatrixFD(:,:,1:nbins))),...
            'xlogaxis',1, 'xlim',[20,20000],'xlabel', 'Frequency (Hz)', ...
            'ylabel', 'Magnitude (dB)', 'stemFlag', 0, 'colors', colors(k,:) ...
            ,'save', save_flag, 'ylim', [-50,5], 'plotYLabels', 0); grid on;
     else
       [plotAxes, plotHandles] = plotImpulseResponseMatrixStem((w/pi)*(fs/2), ...
           20*log10(abs(couplingMatrixFD(:,:,1:nbins))), 'xlogaxis',1, 'xlim',[20,20000],...
           'xlabel', 'Frequency (Hz)', 'ylabel', 'Magnitude (dB)', 'stemFlag', 0, ...
           'colors', colors(k,:) ,'save', save_flag, 'ylim', [-50,5], 'plotYLabels', 1); grid on;
     end
   
     lgdstr{k} = ['a = ', num2str(round(aperture(k),3)), ' m'];
end
set(plotAxes,'XScale','log');


Lgnd = legend(plotAxes(1),lgdstr);
Lgnd.NumColumns = ceil(length(lgdstr)/2);
Lgnd.Position(1:2) = [0.10 0.48];


end