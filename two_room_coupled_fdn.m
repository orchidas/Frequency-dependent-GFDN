%% GFDN with two coupled rooms and filter feedback matrix
%PRODUCES FIGURE 4 AND 5 IN PAPER

close all; clc;
addpath('utils/.');
addpath('plot/.');
addpath(genpath('fdnToolbox/.'));


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
subplot(121);
w = linspace(0,pi,1000);
w = w(1:end-1);
[b1,a1] = shelfeq(wt1,[r1t60_0 r1t60_pi]);
[b2,a2] = shelfeq(wt2,[r2t60_0 r2t60_pi]);
H1 = freqz(b1,a1,w);
H2 = freqz(b2,a2,w);
loglog((w/pi) * (fs/2),abs(H1), 'LineWidth', 1.2, 'LineStyle','-.');hold on;grid on;
loglog((w/pi) * (fs/2),abs(H2), 'LineWidth', 1.2);hold off;
ylabel('T60 (s)');
yticks([0 0.1 0.5 1 2 3 3.5]);
xlabel('Frequency (Hz)');
xlim([20 30000]);
set(gca, 'FontUnits','points', 'FontWeight','normal', 'FontSize',8, 'FontName','Times');


%% FDN attenuation filters

a = zeros(ndel,2);
b = zeros(ndel,2);
col = [[0,0.447,0.741];[0.850,0.325,0.098];[0.929,0.694,0.125];[0 0 0]];
plt = gobjects(1,ndel);
ls = {'-.','-'};

for i = 1:ndel
    gamma_0 = exp(log(0.001)*(tau(i)/fs)/rt60_0(ceil(i/nset)));
    gamma_pi = exp(log(0.001)*(tau(i)/fs)/rt60_pi(ceil(i/nset)));
    [b(i,:),a(i,:)] = shelfeq(wt(ceil(i/nset)),[gamma_0;gamma_pi]);
    H = freqz(b(i,:),a(i,:),w);  
    subplot(122);
    plt(i) = semilogx((w/pi) * (fs/2),20*log10(abs(H)), 'Color',col(ceil(i/nset),:), ...
        'LineWidth', 1.2, 'LineStyle',ls{ceil(i/nset)});
    grid on;hold on;
end
 
ylabel('T60 Filter gain (dB)');
xlabel('Frequency (Hz)');
legend([plt(1), plt(nset+1)],{'Room 1 (small)','Room 2 (big)'},'Location','southwest','FontSize',5);
xlim([20 30000]);ylim([-8 0]);
set(gca, 'FontUnits','points', 'FontWeight','normal', 'FontSize',8, 'FontName','Times');
% print(['figures/coupled_T60_filters_',des,'.eps'], '-depsc');


%% individual room mixing matrices and IR

theta = [1 1]*pi/4;
type = 'none';

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
% audiowrite('audio/r1.wav',h_r1./max(abs(h_r1)),fs);
% audiowrite('audio/r2.wav',h_r2./max(abs(h_r2)),fs);

%% loop through different aperture sizes

aperture = 0.25:0.05:0.5; 
room_dims = [[4,2.5,3];[6,5.5,4.5]];    %dimensions of the coupled rooms
[absorp_coef, area]  = calculate_absorption_coefficient(nGrp, room_dims, rt60_0);

% surface_areas = [5,10];      %areas of two rooms
c = 343;    %speed of sound in m/s
type = 'diffract';  %mixing matrix type
des = 'exact';  %design criteria for diffraction filter

plot_ir = 0;
plot_2_stage = 1;


%% PRODUCES FIGURE 4
% plot coupling matrix response
plot_diffraction_matrix_response(2,fs,c,aperture,'areas',area,'absorption_coeffs', absorp_coef);

if plot_ir

    logtmin = 50;
    tscale = 1000;
    time = 0:1/fs:L-1/fs;
    time_off = 0:1/fs:(L+0.01)-1/fs;
    nzeros = length(time_off) - length(time);
    offset = 1.5*(0:length(aperture)-1);
    wtaps = 0.05*fs; %length of NED window (in samples)
    alpha = ones(nsamp,1) * [5,5,5];  %hyperbolic tangent factor
    ned_start = min(tau)+1; %NED is offset by window length and the time of the first reflection
    ned_end = 2*fs;
    
    NameArray = {'Color','LineStyle', 'LineWidth'};
    ValueArray = { [1,0,0],':',0.75; [0,0,0.8],'-',1.0; [.90,.75,0,1], '-.',0.75;};
    fig = figure('Units','inches', 'Position',[0 0 3.29 4.5],'PaperPositionMode','auto');

%% PRODUCES FIGURE 5
elseif plot_2_stage
    
    T60 = zeros(length(aperture),2,2);
    decayRatio = zeros(length(aperture),1);
    leave_out_first = round(0.05*fs);  %leave out first few ms before calculating envelope
    smbeta = 50;  %smoothing window for energy envelope calculation (in ms)
    fig = figure('Units','inches', 'Position',[0 0 6.6 2.9],'PaperPositionMode','auto');
    offset = 0;
end

%% run loop

for i = 1:length(aperture)
    
   %% design diffraction filter

    [b_diff, a_diff, win_len] = design_diffraction_filter(aperture(i),fs,c,des);

    coupling_areas = absorp_coef.*area + (pi*aperture(i)^2).*(ones(nGrp,1)-absorp_coef);

    %% coupled FDN with FDNTB 

    % with filters in the mixing matrix
    couplingMatrixFilter = two_room_coupling_matrix(aperture(i), b_diff, coupling_areas);
    [zAbsorption, feedbackMatrix, inputGain, outputGain, directGain] = ...
        convert_to_FDNTB_params(nGrp, nSize, b, a, theta, b_drive, c_drive, d, ...
        'couplingMatrix', couplingMatrixFilter);
    h_coup_freq = dss2impz(nsamp, tau.', feedbackMatrix, inputGain, outputGain,...
        directGain,'absorptionFilters', zAbsorption);
    h_coup_freq = h_coup_freq(:,lis_pos);

    
     % no filters in the mixing matrix
    couplingMatrixScalar = reshape([1, (pi*aperture(i)^2)/coupling_areas(1); ...
        -(pi*aperture(i)^2)/coupling_areas(2), 1], [2,2,1]);
    [zAbsorption, feedbackMatrix, inputGain, outputGain, directGain] = ...
        convert_to_FDNTB_params(nGrp, nSize, b, a, theta, b_drive, c_drive, d, ...
        'couplingMatrix', couplingMatrixScalar);
    h_coup = dss2impz(nsamp, tau.', feedbackMatrix, inputGain, outputGain,...
        directGain,'absorptionFilters', zAbsorption);
    h_coup = h_coup(:,lis_pos);

   
    if plot_ir
        %% plot IR
         
        % h_coup is too small, so to make it visible, I am scaling it by
        % the norm of h_r1
        ys = [circshift(h_r2, tau(1)), circshift(h_r1, tau(nset+1)), ...
            circshift(h_coup_freq.*(norm(h_r2)/norm(h_coup_freq)), tau(1)-tau(nset+1))];
        [ned_freq, rep] = rrecho(h_coup_freq, wtaps);
        [ned, rep] = rrecho(h_coup, wtaps);
        figure(fig);

        hfig = semilogx(tscale*(time_off+logtmin/1000), [zeros(nzeros,3);tanh(alpha.*ys)./tanh(alpha)] ...
            + offset(i)); grid on; hold on;
        set(hfig, NameArray, ValueArray);

        %plot NEDs with and without filter feedback matrix
        semilogx(tscale*(time_off(ned_start:ned_end)+logtmin/1000), ...
            circshift(ned_freq(ned_start:ned_end),tau(1))+offset(i),'k', 'LineWidth', 0.5);hold on;
        semilogx(tscale*(time_off(ned_start:ned_end)+logtmin/1000), ...
            circshift(ned(ned_start:ned_end),tau(1))+offset(i),'k', 'LineWidth', 0.75,...
            'LineStyle', '--');hold on;

        text(400,1.5*i-1, ['aperture = ' num2str(round(aperture(i),3)),' m'], 'FontSize',7);
        ylabel('Amplitude');
        xlabel('Time (ms)');
        ylim([-1,1.5*(length(aperture)+1)]);
        xlim(tscale*[logtmin/1000+0.02 logtmin/1000+L-1]);
        legend('Larger room', 'Smaller room', 'Coupled rooms','NED with FFM', ...
            'NED with SFM', 'Location','NW','FontSize',5);
        drawnow;

    %% PRODUCES FIGURE 5
    % fit and plot 2-stage decay
    elseif plot_2_stage
        
        h = [h_coup, h_coup_freq];
        for j = 1:2
            
            ir = h(leave_out_first:end, j); %leave out first 50 ms;
            time = 0:1/fs:(length(ir)-1)/fs;
            ir_env = energy_envelope(ir,fs,smbeta);
            t_env = linspace(0,length(ir)/fs,length(ir_env)).';
            [taus, gamma, ir_fit] = fit_two_stage_decay(ir_env, t_env,[rt60_pi, rt60_0]);
            T60(i,:,j) = taus * log(1000);
            c1 = log(ir_fit(1));     %calculate the right y-intercept assuming slope is right
            line1 = c1 - t_env./taus(1) ;
            c2 = log(ir_fit(end)) + t_env(end)/taus(2);
            line2 = c2 - t_env./taus(2) ;
            tp_x = (c2-c1)/(1/taus(2) - 1/taus(1));
            tp_y = c1 - tp_x/taus(1);
            
            
            %% plot 2 stage decay with and without FFM
            figure(j+101);
%             subplot(1,2,j)
            %if turning point is negative, single slope decay is observed
            if tp_x < 0
                hfig = plot(t_env, ([mag2db(ir_env) mag2db(ir_fit) ...
                    20*log10(exp(1))*line2]+offset));
                hold on;grid on;
                set(hfig,{'Linestyle','Color','LineWidth'},{'-',col(2,:),1.3;...
                    '--',col(3,:),1.3;':',col(4,:),1});

            else
                hfig = plot(t_env, ([mag2db(ir_env) mag2db(ir_fit) ...
                    20*log10(exp(1))*line1 20*log10(exp(1))*line2]+offset));
                hold on;grid on;
                set(hfig,{'Linestyle','Color','LineWidth'},{'-',col(2,:),1.3;...
                    '--',col(3,:),1.3;':',col(4,:),1;':',col(4,:),1});

            end
            plot(tp_x, 20*log10(exp(1))*tp_y + offset, 'bo', 'MarkerSize',5, 'MarkerFaceColor','b');
            hold on;
            xpos = find(t_env >= 0.5, 1, 'first');
            text(t_env(xpos),mag2db(ir_env(xpos))+offset+10, ['a =' num2str(round(aperture(i),3)),' m'],...
                'FontSize',7, 'interpreter','latex');
            drawnow;
        end
        offset = offset+40;

    end
    %% play and save IR
%      audiowrite(['audio/coupled_aperture=',num2str(round(aperture(i),3)),...
%                 '_freq_', des,'.wav'], h_coup_freq./max(abs(h_coup_freq)),fs);
%      audiowrite(['audio/coupled_aperture=',num2str(round(aperture(i),3)),...
%                  '_', des,'.wav'], h_coup./max(abs(h_coup)),fs);
%      soundsc(h_coup,fs); pause(2); soundsc(h_coup_freq,fs);

end

%% save figures 

if plot_ir
    % coupled RIR
    off = [offset; offset+1];
    yticks(off(:)'); yticklabels({'0','1','0','1','0','1','0','1','0','1'});
    set(gca, 'FontUnits','points', 'FontWeight','normal', 'FontSize',8, 'FontName','Times');
%     print(['figures/coupled_FDN_IR_',des,'.eps'], '-depsc');

%% PRODUCES FIGURE 5
elseif plot_2_stage
    % two stage decay

%     ttl = {'Scalar feedback matrix', 'Filter feedback matrix'};
    for j = 1:2
%         subplot(1,2,j);
        figure(j+101)
        set(gcf,'Units','inches', 'Position',[0 0 2.29 4.5],'PaperPositionMode','auto');
        hold off;
        ylim([-120, 180]);xlim([-0 2]); 
        xlabel('Time (s)');
        ylabel('Magnitude (dB)');
%         title(ttl{j});
        set(gca, 'FontUnits','points', 'FontWeight','normal', 'FontSize',8, 'FontName','Times');
        saveas(gcf,['figures/coupled_2_stage_decay_',num2str(j),des,'.png']);
    end
%     print(['figures/coupled_2_stage_decay_',des,'.eps'], '-depsc');


end


