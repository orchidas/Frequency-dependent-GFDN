%% Multi-room coupling example with Frequency-dependent GFDNs
% PRODUCES FIGURE 8 IN PAPER

close all; 
addpath('../.');
addpath('../plot/.')
addpath('../utils/.');
addpath(genpath('../fdnToolbox/.'));

%sampling frequency
fs = 48000;
%number of delay lines = 12
nDel = 12;
%input gains - speaker in first room
b_drive = [-0.7429   -0.8463   -0.2649   -0.8975   zeros(1,8)]';
%output gains - mic in third room
c_drive = [zeros(8,1); ones(4,1)];
%number of groups in GFDN
nGrp = 3;   
%direct path gain
d = 0;
%number of delay lines for each room
nSize = [4,4,4];   
assert(sum(nSize) == nDel);
% whether to save RIR
saveRIR = false;

%% delay line lengths (co-prime)
tau = sort([1023, 511, 897, 1559, 617, 987, 781, 2029, 613, 1443, 2393, 1777])';

%% T60 filters for coupled rooms

% DC T60
rt60_0 = [0.8;1.5;3.5];
% Nyquist T60
rt60_pi = [0.2;0.8;1.8];
% transition frequency
wt = [1000;2000;4000]/(fs/2);


%% PRODUCES FIGURE 8A 
% Desired T60 response of nGrp rooms

fig = figure('Units','inches', 'Position',[0 0 3.29 2.2],'PaperPositionMode','auto');

[b1,a1] = shelfeq(wt(1),[rt60_0(1) rt60_pi(1)]);
[b2,a2] = shelfeq(wt(2),[rt60_0(2) rt60_pi(2)]);
[b3,a3] = shelfeq(wt(3),[rt60_0(3) rt60_pi(3)]);

plot_desired_T60(fs, nGrp, [b1; b2; b3], [a1; a2; a3], fig);
% exportgraphics(gcf, '../figures/t60_rooms.pdf','BackgroundColor','none','ContentType','vector');

%% FDN attenuation filters

a = zeros(nDel,2);
b = zeros(nDel,2);
col = [[0,0.447,0.741];[0.850,0.325,0.098];[0.929,0.694,0.125];[0 0 0]];
ls = {'-.','-',':'};

fig = figure('Units','inches', 'Position',[0 0 3.29 2.2],'PaperPositionMode','auto');

for i = 1:nDel
    gamma_0 = exp(log(0.001)*(tau(i)/fs)/rt60_0(ceil(i/nSize(1))));
    gamma_pi = exp(log(0.001)*(tau(i)/fs)/rt60_pi(ceil(i/nSize(1))));
    [b(i,:),a(i,:)] = shelfeq(wt(ceil(i/nSize(1))),[gamma_0;gamma_pi]);
end
 
plot_fdn_attenuation_filters(nDel, nSize, fs, b, a, col, ls);
% exportgraphics(gcf, '../figures/multiroom_t60_filters.pdf','BackgroundColor','none','ContentType','vector');

%% individual room mixing matrices and IR

% mixing angle for individual rooms (Hadamard)
theta = [1 1 1]*pi/4;

%individual room FDN IRs
L = 2;   %length of impulse response to be calculated
nsamp = L*fs;

h_rooms = zeros(nsamp, nGrp);
start = 0;
for i = 1:nGrp

    [zAbsorption, feedbackMatrix, inputGain, outputGain, directGain] = ...
        convert_to_FDNTB_params(1, nSize(i), ...
        b(start+1:start+nSize(i),:), a(start+1:start+nSize(i),:), ...
        theta(i), b_drive(1:nSize(i)), ones(nSize(i),1), d);
    
    h_rooms(:,i) = dss2impz(nsamp, tau(start+1:start+nSize(i)).', feedbackMatrix, ...
        inputGain, outputGain, directGain,'absorptionFilters', zAbsorption);
    
    start = start + nSize(i);

    % save results
    if saveRIR
        audiowrite(['../audio/multiroom/r', num2str(i),'.wav'],h_rooms(:,i)./max(abs(h_rooms(:,i))),fs);
    end
end

%% coupled rooms IR

% aperture size connecting the three rooms - assumed to be symmetric
aperture = [0, 0.22, 0.35; 0.22 0 0.5; 0.35 0.5 0];
% room dimensions
room_dims = [[3,1.8,2]; [4,3.2,3.8]; [6,5.5,4.5]];
% calculate absorption coefficients from T60 and room dimensions
[absorp_coef, area]  = calculate_absorption_coefficient(nGrp, room_dims, rt60_0);
des = 'exact';  %design criteria for diffraction filter
c = 343;        %speed of sound in m/s

%use FDNTB
absorption.b = reshape(b, [nDel,1,2]);
absorption.a = reshape(a, [nDel,1,2]);
zAbsorption = zTF(absorption.b,absorption.a,'isDiagonal',true);

% block mixing matrix
M_block = get_block_mixing_matrix(nSize, nGrp, theta);

% filter coupling matrix
[couplingMatrix, ~, ~] = multiroom_coupling_matrix(fs, c, nGrp, aperture, absorp_coef, area, des);

%% PRODUCES FIGURE 8B
% plot T60 deviation

fig = figure('Units','inches', 'Position',[0 0 3.29 2.2],'PaperPositionMode','auto');
plot_T60_deviation(fs, couplingMatrix, tau, [b1; b2; b3], [a1; a2; a3], col, fig); 
exportgraphics(gcf, '../figures/t60_error_bounds.pdf','BackgroundColor','none','ContentType','vector');

%% generate impulse response

% feedback matrix
degree = max(max(win_len));
feedbackMatrix = zeros(nDel, nDel, degree);
for k = 1:degree
    feedbackMatrix(:,:,k) =  kron(couplingMatrix(:,:,k), ones(nDel/nGrp, nDel/nGrp)) .* M_block;
end

% find impulse response
outputGain = blkdiag(ones(1,nSize(1)), ones(1,nSize(1)), ones(1,nSize(1)));
direct = zeros(nGrp,1);
h_coup_freq_fdntb = dss2impz(nsamp, tau.', real(feedbackMatrix), b_drive, outputGain, direct,...
    'absorptionFilters', zAbsorption);

%% PRODUCES FIGURE 8C
% plot spectrograms

ftgram2N4S([h_rooms, h_coup_freq_fdntb(:,3)./max(abs(h_coup_freq_fdntb(:,3)))], ...
    fs, '../figures/multiroom_spec_causal', true);


%% save impulse responses

if saveRIR
    audiowrite('../audio/multiroom/s=1_m=1_causal_fdntb.wav',...
                h_coup_freq_fdntb(:,1)./max(abs(h_coup_freq_fdntb(:,1))),fs);
    audiowrite('../audio/multiroom/s=1_m=2_causal_fdntb.wav',...
                h_coup_freq_fdntb(:,2)./max(abs(h_coup_freq_fdntb(:,2))),fs);
    audiowrite('../audio/multiroom/s=1_m=3_causal_fdntb.wav',...
                h_coup_freq_fdntb(:,3)./max(abs(h_coup_freq_fdntb(:,3))),fs);
end

