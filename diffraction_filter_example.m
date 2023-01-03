%% Design diffraction filter in the frequency domain and invert it
% PRODUCES FIGURE 3 IN PAPER
close all;
clear all;

%% physical parameters and diffraction model

% aperture = 0.05:0.25:1.25;   % size/radius of aperture, 10cm-2m
aperture = 5*logspace(-2,-1,5); %radius of aperture
% aperture = 5*logspace(-3,0,6); %radius of aperture
c = 343;    %speed of sound
Fs = 44100;    %sampling rate
ntheta = 50;
theta = linspace(0, pi/2 - 0.001, ntheta);


% arg = aperture(i)/pi * (k'*sin(theta));
% tau = sigma * ((sinc(arg).^2)./(ones(nfreq,1) * cos(theta)));
% maxTau = max(tau(:,end));
% 
% figure;
% colormap(jet);
% [x,y] = meshgrid(theta, freqs);
% h = surf(x,y,tau);
% set(h, 'LineStyle', 'none');
% zlabel('Transmission coefficient, \tau');
% xlabel('Incidence angle, \theta');
% ylabel('Frequency in Hz');


%% loop over different aperture sizes

fig = figure('Units','inches', 'Position',[0 0 3.3 4.5],'PaperPositionMode','auto');
% fig = figure(1);
col = get(gca,'colororder');
NameArray = {'LineStyle'};
ValueArray = {'-';'--'};
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
    
    
    sigma = pi * aperture(i)^2;  %assuming circular aperture
    fun = @(angle) sigma * (sinc((k * aperture(i))/pi * sin(angle)).^2).*sin(angle);
    % fun = @(angle) (sigma/cos(angle)) * sinc((k * aperture(i))/pi * sin(angle)).^2;
    tau_avg = 2*integral(fun, 0 , pi/2, 'ArrayValued',true); 

    % we know that the transmission coefficient should be 1 for high
    % frequencies
    cutoff = find(k*aperture(i) >= 5,1,'first');
    tau_avg(cutoff:end) = 1.0;
    startFrom = find(k*aperture(i) >= 3, 1, 'first');
    %interpolation to avoid abrupt jump
    tau_avg(startFrom:cutoff) = linspace(tau_avg(startFrom), tau_avg(cutoff), length(startFrom:cutoff));

%   use Parks mcClellan to design FIR filter
%     F = [0; freqs(cutoff-1)/(Fs/2); freqs(cutoff+1)/(Fs/2); 1];
%     A = [tau_avg(1); tau_avg(cutoff); tau_avg(end)];
%     [win_len, f0, a0, w] = firpmord(F.',A.',[0.1, 0.01, 0.1],fs);
    tau_pm = firpm(win_len-1, freqs/(Fs/2), tau_avg); 
    [H_pm,~] = freqz(tau_pm,1,nfreq);

    % inverse FFT from integral - this is what filter should be
    tau_spec = [tau_avg, conj(tau_avg(2:end-1))];
    % a real signal will have a hermitian symmetric spectrum
    tau_filt = fftshift(ifft(tau_spec, 'symmetric'));
    %area of tau_filt remains same, regardless of aperture size
    area = 0.5 * win_len * tau_filt(nfreq);

    a = 1;
    b = win;
    % normalize the window to have proper magnitude
    b = b .* (Fs*sigma/(4*c)*pi/win_len);
    [H,w] = freqz(b,a,nfreq);
    
    %% PRODUCES FIGURE 3
    fig;
    ax1 = subplot(211);
    h1 = semilogx(freqs, 20*log10(abs(tau_avg)), 'LineStyle','-');grid on;hold on;
    h2 = semilogx(freqs, 20*log10(abs(H_pm)), 'LineStyle','-');grid on;hold on;
    set([h1,h2], 'Color', col(i,:), NameArray, ValueArray);
    xlabel('Frequency (Hz)');xlim([20, Fs/2]);
    ylabel('$\zeta_{avg}(f)$ magnitude (dB)', 'Interpreter','latex');
   
    ax2 = subplot(212);
    time_aligned_tau = tau_filt(nfreq-(win_len-1)/2 : nfreq+(win_len-1)/2); 
    h = plot([time_aligned_tau.', tau_pm.'], 'Color', col(i,:));grid on; hold on;
    axis tight;
    xlabel('Number of samples, n'); ylabel('$\zeta_{avg}(n)$ amplitude', 'Interpreter','latex');
    lgdstr{i} =  ['a = ', num2str(round(aperture(i),3)),' m'];
    leglines = [leglines, h(1)];
    set(h, NameArray, ValueArray);
    
end

fig;hold off;
set([ax1,ax2],'FontUnits','points', 'FontWeight','normal', 'FontSize',8, 'FontName','Times');
legend(leglines,lgdstr);
print('figures/diffraction_filter.eps', '-depsc');



