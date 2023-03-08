function [] = plot_coupled_rir(fs, L, tau, aperture, h_r1, h_r2, h_coup_scalar, h_coup_filter, offset, i, nset, tanh_factor)
%% 
% Plot coupled RIR and NED for varying aperture lengths
% Inputs:
% fs - sampling rate
% L - length of FDN IR in seconds
% tau - delay line lengths in samples
% aperture - aperture size
% h_r1, h_r2 - impulse responses of decoupled rooms
% h_coup_scalar - impulse response with scalar feedback matrix
% h_coup_filter - impulse response with filter feedback matrix
% i - counter variable value (if evaluating for many aperture sizes)
% offset - by how much is each plot offset
% nset - number of delay lines in group 1
% tanh_factor - hyperbolic tangent factor for normalisation


    nsamp = L*fs;
    logtmin = 50;
    tscale = 1000;
    time = 0:1/fs:L-1/fs;
    time_off = 0:1/fs:(L+0.01)-1/fs;
    nzeros = length(time_off) - length(time);
    wtaps = 0.02*fs; %length of NED window (in samples)
    alpha = ones(nsamp,1) * (ones(1,3) * tanh_factor);  %hyperbolic tangent factor
    ned_start = min(tau)+1; %NED is offset by window length and the time of the first reflection
    ned_end = nsamp;

    NameArray = {'Color','LineStyle', 'LineWidth'};
    ValueArray = { [1,0,0],':',0.75; [0,0,0.8],'-',1.0; [.90,.75,0,1], '-.',0.75;};

    % h_coup is too small, so to make it visible, I am scaling it by
    % the norm of h_r1
    ys = [circshift(h_r2, tau(1)), circshift(h_r1, tau(nset+1)), ...
        circshift(h_coup_filter.*(norm(h_r2)/norm(h_coup_filter)), tau(1)-tau(nset+1))];
    [ned_freq, ~] = rrecho(h_coup_filter, wtaps);
    [ned_scalar, ~] = rrecho(h_coup_scalar, wtaps);

    % plot RIRs
    figure(101);
        set(gcf,'Units','inches', 'Position',[0 0 2.29 4.5],'PaperPositionMode','auto'););
    hfig = semilogx(tscale*(time_off+logtmin/1000), [zeros(nzeros,3);tanh(alpha.*ys)./tanh(alpha)] ...
        + offset(i)); grid on; hold on;
    set(hfig, NameArray, ValueArray);

    %plot NEDs with scalar and with filter feedback matrix
    semilogx(tscale*(time_off(ned_start:ned_end)+logtmin/1000), ...
        circshift(ned_freq(ned_start:ned_end),tau(1))+offset(i),'k', 'LineWidth', 0.5);hold on;
    semilogx(tscale*(time_off(ned_start:ned_end)+logtmin/1000), ...
        circshift(ned_scalar(ned_start:ned_end),tau(1))+offset(i),'k', 'LineWidth', 0.75,...
        'LineStyle', '--');hold on;

    text(400,1.5*i-1, ['aperture = ' num2str(round(aperture(i),3)),' m'], 'FontSize',7);
    ylabel('Amplitude');
    xlabel('Time (ms)');
    ylim([-1,1.5*(length(aperture)+1)]);
    xlim(tscale*[logtmin/1000+0.02 logtmin/1000+L-1]);
    legend('Larger room', 'Smaller room', 'Coupled rooms','NED with FFM', ...
        'NED with SFM', 'Location','NW','FontSize',5);
    drawnow;

end

