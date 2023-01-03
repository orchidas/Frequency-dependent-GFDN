%% Paraunitary coupling matrix in 3 coupled rooms
% PRODUCES FIGURES 6 AND 7 IN PAPER
addpath('plot/.');
addpath('utils/.');
close all;

%% design coupling matrix

Nrooms = 3;
% aperture size connecting the three rooms - assumed to be symmetric
% aperture = [0, 0.08, 0.05; 0.08 0 0.1; 0.05 0.1 0];
aperture = [0, 0.22, 0.35; 0.22 0 0.5; 0.35 0.5 0];
room_dims = [[3,1.8,2]; [4,3.2,3.8]; [6,5.5,4.5]];
rt60_0 = [0.8;1.5;3.5];
rt60_pi = [0.2;0.8;1.8];
[absorp_coef, area]  = calculate_absorption_coefficient(Nrooms, room_dims, rt60_0);
coupling_areas = ones(Nrooms, Nrooms);

fs = 44100;
c = 343;    %speed of sound in m/s
type = 'diffract';  %mixing matrix type
des = 'pm';  %design criteria for diffraction filter
filter_coefs = zeros(Nrooms,Nrooms,300);
win_len = zeros(Nrooms,Nrooms);
Nfft = 1024;
w = linspace(0, pi, Nfft/2);


for i = 1:Nrooms

    for j = 1:Nrooms
       coupling_areas(i,j) = absorp_coef(i).*area(i) + (pi*aperture(i,j)^2).*(1-absorp_coef(i));

        if (i == j)
            continue;
        else
            [b_diff, a_diff, win_len(i,j)] = design_diffraction_filter(aperture(i,j),fs,c,des);
            filter_coefs(i,j,1:win_len(i,j)) =  b_diff;

        end
    end
end
degree = max(max(win_len));
couplingMatrix = zeros(Nrooms, Nrooms, degree);
const = max(max(max(filter_coefs)));
zeroPos = 1; %(degree+1)/2;     %position of 0th order coefficient
midPoint = (degree+1)/2;
signVec = zeros(1,1,degree);


for i = 1:Nrooms
    for j = 1:Nrooms
        if ( i == j)
            continue;
        else
            shift = (win_len(i,j)+1)/2;
            couplingMatrix(i,j,midPoint-shift+1:midPoint+shift-1) = (pi*(aperture(i,j)^2))/sqrt(coupling_areas(i,j)*coupling_areas(j,i)) ...
                .* filter_coefs(i,j,1:win_len(i,j));
        end
    end
end


%% find closest PU matrix with per frequency Procrustes

[couplingMatrix, originalMatrix, signComb] = make_more_paraunitary(couplingMatrix);
couplingMatrix_freq = fft(couplingMatrix, Nfft,3);
originalMatrix_freq = fft(originalMatrix, Nfft, 3);



%% PRODUCES FIGURE 6A
% plot time domain filters
figure('Units','inches', 'Position',[0 0 2.9 2.5],'PaperPositionMode','auto');
[~, pHDes] = plotImpulseResponseMatrix(0:(degree-1), originalMatrix, 'save', 1, 'colors', [0, 0.4470, 0.7410]...
   ,'ylim', [-0.01,0.01], 'title', 'Initial matrix', 'stemFlag', 0, 'xlabel', 'Samples', 'ylabel', 'Amplitude'); hold on;
print('figures/3x3_couplingMatrix_init.eps', '-depsc');


figure('Units','inches', 'Position',[0 0 2.9 2.5],'PaperPositionMode','auto');
[~, perFreqProcrustes] = plotImpulseResponseMatrix(0:(degree-1), couplingMatrix, 'save', 1, 'colors', [0.8500, 0.3250, 0.0980]...
   ,'ylim', [-0.01,0.01] , 'title', 'Per-Frequency Procrustes', 'stemFlag',0, 'xlabel', 'Samples', 'ylabel', 'Amplitude'); 
print('figures/3x3_couplingMatrix_pFProc.eps', '-depsc');


%% PRODUCES FIGURE 7
% plot frequency response

figure('Units','inches', 'Position',[0 0 6.6 2.5],'PaperPositionMode','auto');
[pADes, pHDes] = plotImpulseResponseMatrix((w/pi)*(fs/2), 20*log10(abs(originalMatrix_freq(:,:,1:Nfft/2))), 'xlogaxis',1, 'xlim',[20,20000],...
    'xlabel', 'Frequency (Hz)', 'ylabel', 'Magnitude (dB)', 'stemFlag',0, 'save', 1, 'plotYLabels', 0); hold on;
[pAProcrustes, pProcrustes] = plotImpulseResponseMatrix((w/pi)*(fs/2), 20*log10(abs(couplingMatrix_freq(:,:,1:Nfft/2))), 'xlogaxis',1, 'xlim',[20,20000],...
     'xlabel', 'Frequency (Hz)', 'ylabel', 'Magnitude (dB)', 'stemFlag',0, 'save',1, 'commonYAxis', 1, 'ylim', [-80,5]); hold off;


% add legend
Lgnd = legend(pADes(1,1), [pHDes(1,1), pProcrustes(1,1)], 'Desired', 'PF Procrustes', 'Location', 'southeast');
print('figures/3x3_couplingMatrix_resp_compared.eps', '-depsc');


%% PRODUCES FIGURE 6B
% plot magnitude of singular values (SVD per frequency bin)

color_gradient = zeros(Nrooms, 3, 2);
blue = [0 0 1];
light_blue = [91, 207, 244] / 255;
color_gradient(:,:,1) = create_color_gradient(blue, light_blue, Nrooms);

red = [204 0 0]/255;
pink = [255, 192, 203]/255;
color_gradient(:,:,2) = create_color_gradient(red, pink, Nrooms);


figure('Units','inches', 'Position',[0 0 2.9 2.5],'PaperPositionMode','auto'); grid on;
S = svdPerBin(fft(originalMatrix,2*degree,3));
p1 = semilogx(linspace(-fs/2,fs/2,2*degree), 20*log10(abs(S')));grid on;
xlabel('Frequency (Hz)');
ylabel('Singular Values (dB)');
NameArray = {'Color', 'LineStyle'};
ValueArray = {color_gradient(1,:,1), '-'; color_gradient(2,:,1), '--'; color_gradient(3,:,1), '-.'};
set(p1, NameArray, ValueArray);
xlim([0,fs/2]);
axis tight;
set(gca,'FontUnits','points', 'FontWeight','normal', 'FontSize',8, 'FontName','Times');
print('figures/3x3_couplingMatrix_SVD_init.eps', '-depsc');


figure('Units','inches', 'Position',[0 0 2.9 2.5],'PaperPositionMode','auto'); grid on;
S = svdPerBin(fft(couplingMatrix,2*degree,3));
p2 = semilogx(linspace(-fs/2,fs/2,2*degree), 20*log10(abs(S'))); grid on;
xlabel('Frequency (Hz)'); ylabel('Singular Values (dB)');
NameArray = {'Color', 'LineStyle'};
ValueArray = {color_gradient(1,:,2), '-'; color_gradient(2,:,2), '--'; color_gradient(3,:,2), '-.'};
set(p2, NameArray, ValueArray);
axis tight;
xlim([0, fs/2]);
set(gca,'FontUnits','points', 'FontWeight','normal', 'FontSize',8, 'FontName','Times');
print('figures/3x3_couplingMatrix_SVD_pFProc.eps', '-depsc');


%% helper functions

function PP = svdPerBin(TT)

    PP = zeros(size(TT,[1 3]));
    for f = 1:size(TT,3)
       s = (svd(TT(:,:,f)));  
       PP(:,f) = s;
    end
end


function [color_gradient] = create_color_gradient(col1, col2, nterms)
    color_gradient = [linspace(col1(1),col2(1),nterms)', linspace(col1(2),col2(2),nterms)'...
        , linspace(col1(3),col2(3),nterms)'];
end
