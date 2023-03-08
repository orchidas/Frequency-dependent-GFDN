%% Design diffraction filter in the frequency domain and invert it
% PRODUCES FIGURE 3 IN PAPER
close all;
clear all;

%% physical parameters and diffraction model

aperture = 5*logspace(-2,-1,5); %radius of aperture
c = 343;    %speed of sound
Fs = 44100;    %sampling rate
ntheta = 50;
theta = linspace(0, pi/2 - 0.001, ntheta);

%% loop over different aperture sizes

fig1 = figure('Units','inches', 'Position',[0 0 3.3 2.4],'PaperPositionMode','auto'); grid on; hold on;
fig2 = figure('Units','inches', 'Position',[0 0 3.3 2.0],'PaperPositionMode','auto');
col = get(gca,'colororder');
NameArray = {'LineStyle'};
ValueArray = {'-';'-.'};
lgdstr = {};
leglines = [];

for i = 1:length(aperture)
    
    %% %% design FIR filter to match diffraction response
    order = 4*(aperture(i) * Fs/c);
    win_len = round(order);

    % make sure window length is odd
    if (mod(win_len,2) == 0)
        win_len = win_len + 1;
    end

    win = bartlett(win_len);
    nfreq = 2^12;
    freqs = linspace(0,Fs/2, nfreq);
    k = 2*pi*freqs/c;   %wave number


    %% numerically integrate over theta to get average energy arriving at aperture
    
    % use Parks mcClellan to design FIR filter
    tau_pm = design_diffraction_filter(aperture(i), Fs, c, 'pm');
    [H_pm,~] = freqz(tau_pm,1,nfreq);

    % inverse FFT from integral - this is what filter should be
    [tau_filt, a, len, tau_avg] = design_diffraction_filter(aperture(i), Fs, c, 'exact');
%     %area of tau_filt remains same, regardless of aperture size
%     area = 0.5 * win_len * tau_filt(nfreq);

    
    %% PRODUCES FIGURE 3
    set(0, 'currentfigure', fig1);
    h = plot([tau_filt, tau_pm], 'Color', col(i,:)); 
    axis tight;
    xlabel('Time (samples)'); ylabel('$\zeta_{avg}(n)$ amplitude', 'Interpreter','latex');
    lgdstr{i} =  ['a = ', num2str(round(aperture(i),3)),' m'];
    leglines = [leglines, h(1)];
    set(h, NameArray, ValueArray);
    
    set(0, 'currentfigure', fig2);
    h1 = semilogx(freqs, 20*log10(abs(tau_avg)), 'LineStyle','-');grid on;hold on;
    h2 = semilogx(freqs, 20*log10(abs(H_pm)), 'LineStyle','-');grid on;hold on;
    set([h1,h2], 'Color', col(i,:), NameArray, ValueArray);
    xlabel('Frequency (Hz)');xlim([20, Fs/2]);
    ylabel('$\zeta_{avg}(f)$ magnitude (dB)', 'Interpreter','latex');
   
end

set(0, 'currentfigure', fig1); hold off;
set(gca,'Position',[0.1300 0.1100 0.7750 0.7350])
Lgnd = legend(leglines,lgdstr,'Location','north');
Lgnd.NumColumns = ceil(length(lgdstr)/2);
Lgnd.Position(1) = 0.03;
Lgnd.Position(2) = 0.88;
exportgraphics(gcf, '../figures/diffraction_filter1.pdf','BackgroundColor','none','ContentType','vector');


set(0, 'currentfigure', fig2); hold off;
exportgraphics(gcf, '../figures/diffraction_filter2.pdf','BackgroundColor','none','ContentType','vector');

