%% GFDN with two coupled rooms and filter feedback matrix
%PRODUCES FIGURE 4 AND 5 IN PAPER

close all; clc;
addpath('../utils/.');
addpath(genpath('../plot/.'));
addpath(genpath('../fdnToolbox/.'));


%sampling frequency
fs = 48000;
%number of delay lines = 8
n = 4;
ndel = 2^n;
%number of coupled rooms
nGrp = 2;
%number of delay lines for each room
nset = ndel/nGrp;
nSize = nset .* ones(nGrp,1);
%input gains
% source in small room 
b_drive = [-0.9002, 0.9774, 0.2956, -0.6900, -0.2238, 0.2147, 0.4905, -0.6664, zeros(1,nset)].'; 
%output gains
c_drive = [zeros(nset,1); ones(nset,1)];  

% source and listener positions
src_pos = 1;
lis_pos = 2;

%direct path gain
d = 0;

% whether to save RIRs
saveRIR = false;

% what plots to save
plot_ir = false;
plot_2_stage = true;

%% delay line lengths

% first set associated with smaller room, second set associated with bigger
tau = sort([743, 567, 919, 869, 1667, 2023, 1227, 1561, 607, 1009, 479, 1361, 1171, 733, 829, 1123])';


%% T60 filters for coupled rooms, room 1 is smaller (smaller T60) than room 2 

r1t60_0 = 2;
r2t60_0 = 3;
rt60_0 = [r1t60_0; r2t60_0];

r1t60_pi = 0.5;
r2t60_pi = 1;
rt60_pi = [r1t60_pi; r2t60_pi];

wt1 = 1000/(fs/2);
wt2 = 4000/(fs/2);
wt = [wt1; wt2];


%% T60 response of 2 rooms

fig = figure('Units','inches', 'Position',[0 0 3.29 2.2],'PaperPositionMode','auto');

[b1,a1] = shelfeq(wt1,[r1t60_0 r1t60_pi]);
[b2,a2] = shelfeq(wt2,[r2t60_0 r2t60_pi]);

plot_desired_T60(fs, 2, [b1; b2], [a1;a2], fig)

%% FDN attenuation filters

a = zeros(ndel,2);
b = zeros(ndel,2);

for i = 1:ndel
    gamma_0 = exp(log(0.001)*(tau(i)/fs)/rt60_0(ceil(i/nset)));
    gamma_pi = exp(log(0.001)*(tau(i)/fs)/rt60_pi(ceil(i/nset)));
    [b(i,:),a(i,:)] = shelfeq(wt(ceil(i/nset)),[gamma_0;gamma_pi]);   
end

col = [[0,0.447,0.741];[0.850,0.325,0.098];[0.929,0.694,0.125];[0 0 0]];
ls = {'-.','-'};
plot_fdn_attenuation_filters(ndel, nSize, fs, b, a, col, ls);

%% individual room mixing matrices and IR

theta = [1 1]*pi/4;

%individual room FDN IRs
L = 2;   %length of impulse response to be calculated
nsamp = L*fs;

[zAbsorption, feedbackMatrix, inputGain, outputGain, directGain] = ...
    convert_to_FDNTB_params(1, nset, b(1:nset,:), a(1:nset,:), theta(1), ...
    b_drive(1:nset), ones(nset,1), d);
h_r1 = dss2impz(nsamp, tau(1:nset).', feedbackMatrix, inputGain, outputGain,...
    directGain,'absorptionFilters', zAbsorption);

[zAbsorption, ~, inputGain, outputGain, directGain] = ...
    convert_to_FDNTB_params(1, nset, b(nset+1:end,:), a(nset+1:end,:), theta(2),...
    b_drive(1:nset), ones(nset,1), d);

h_r2 = dss2impz(nsamp, tau(nset+1:end).', feedbackMatrix, inputGain, outputGain, ...
    directGain,'absorptionFilters', zAbsorption);


% save results
if saveRIR
    audiowrite('audio/r1.wav',h_r1./max(abs(h_r1)),fs);
    audiowrite('audio/r2.wav',h_r2./max(abs(h_r2)),fs);
end

%% loop through different aperture sizes

aperture = 0.25:0.05:0.5; 
room_dims = [[4,2.5,3];[6,5.5,4.5]];    %dimensions of the coupled rooms
[absorp_coef, area]  = calculate_absorption_coefficient(nGrp, room_dims, rt60_0);

% surface_areas = [5,10];      %areas of two rooms
c = 343;    %speed of sound in m/s
type = 'diffract';  %mixing matrix type
des = 'exact';  %design criteria for diffraction filter


%% PRODUCES FIGURE 4

% plot coupling matrix response
figure('Units','inches', 'Position',[0 0 3.25 3.3],'PaperPositionMode','auto');
plot_diffraction_matrix_response(2,fs,c,aperture,'areas',area,'absorption_coeffs',...
    absorp_coef, 'save_flag', 1);
exportgraphics(gcf,'../figures/coupled_FDN_couplingMatrix.pdf','BackgroundColor','none','ContentType','vector')


%% PRODUCES FIGURE 5
if plot_ir
  tanh_factor  = 3;
  offset = 1.5*(0:length(aperture)-1);
  
elseif plot_2_stage
  offset = 0;
  all_T60 = zeros(length(aperture), 2, 2);
end

%% run loop

for i = 1:length(aperture)  

    %% coupled FDN with FDNTB 
    
    % with filters in the mixing matrix
    couplingMatrixFilter = two_room_coupling_matrix(nGrp, aperture(i), area, absorp_coef, fs, c, des);
    
    [zAbsorption, feedbackMatrix, inputGain, outputGain, directGain] = ...
        convert_to_FDNTB_params(nGrp, nSize, b, a, theta, b_drive, c_drive, d, ...
        'couplingMatrix', couplingMatrixFilter);
    h_coup_filter = dss2impz(nsamp, tau.', feedbackMatrix, inputGain, outputGain,...
        directGain,'absorptionFilters', zAbsorption);
    h_coup_filter = h_coup_filter(:,lis_pos);

    
     % no filters in the mixing matrix
    coupling_areas = absorp_coef.*area + (pi*aperture(i)^2).*(ones(nGrp,1)-absorp_coef);
    couplingMatrixScalar = reshape([1, (pi*aperture(i)^2)/coupling_areas(1); ...
        -(pi*aperture(i)^2)/coupling_areas(2), 1], [2,2,1]);
    
    [zAbsorption, feedbackMatrix, inputGain, outputGain, directGain] = ...
        convert_to_FDNTB_params(nGrp, nSize, b, a, theta, b_drive, c_drive, d, ...
        'couplingMatrix', couplingMatrixScalar);
    h_coup_scalar = dss2impz(nsamp, tau.', feedbackMatrix, inputGain, outputGain,...
        directGain,'absorptionFilters', zAbsorption);
    h_coup_scalar = h_coup_scalar(:,lis_pos);

   
    if plot_ir
        %% plot IR
        plot_coupled_rir(fig, fs, L, tau, aperture, h_r1, h_r2, h_coup_scalar, h_coup_filter, i, nset, tanh_factor);
        
    %% PRODUCES FIGURE 5
    % fit and plot 2-stage decay
    elseif plot_2_stage
        
        h = [h_coup_scalar, h_coup_filter];
        all_T60(i,:,:) = plot_two_stage_decay(fs, L, aperture, h, rt60_pi, rt60_0, col, i, offset);
        offset = offset+40;

    end
    %% play and save IR
    if saveRIR
         audiowrite(['audio/coupled_aperture=',num2str(round(aperture(i),3)),...
                    '_freq_', des,'.wav'], h_coup_freq./max(abs(h_coup_freq)),fs);
         audiowrite(['audio/coupled_aperture=',num2str(round(aperture(i),3)),...
                     '_', des,'.wav'], h_coup./max(abs(h_coup)),fs);
         soundsc(h_coup,fs); pause(2); soundsc(h_coup_freq,fs);
    end

end

%% save figures 

if plot_ir
    % coupled RIR
    figure(101);
    set(gcf, 'Units','inches', 'Position',[0 0 3.29 4.5],'PaperPositionMode','auto');
    off = [offset; offset+1];
    yticks(off(:)'); yticklabels({'0','1','0','1','0','1','0','1','0','1'});
    set(gca, 'FontUnits','points', 'FontWeight','normal', 'FontSize',8, 'FontName','Times');
    print(['../figures/coupled_FDN_IR_',des,'.eps'], '-depsc');

%% PRODUCES FIGURE 5
elseif plot_2_stage
    % two stage decay
    for j = 1:2
        figure(j+101)
        set(gcf,'Units','inches', 'Position',[0 0 2.29 4.5],'PaperPositionMode','auto');
        hold off;
        set(gca, 'FontUnits','points', 'FontWeight','normal', 'FontSize',8, 'FontName','Times');
        print(['../figures/coupled_2_stage_decay_',num2str(j),des,'.eps'], '-depsc');
    end
end


