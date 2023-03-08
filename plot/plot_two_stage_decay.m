function [T60] = plot_two_stage_decay(fs, L, aperture, h, rt60_pi, rt60_0, col, i, offset)
%% Plot two stage decay of coupled FDN RIRs
% Inputs:
% fs - sampling frequency
% L - length of IR in seconds
% aperture - aperture size
% h - RIR to be analysed
% rt60_pi - limits of Nyquist T60
% rt60_0 - limits of DC T60
% col - colors to be used
% i - counter variable value (if evaluating for many aperture sizes)
% offset - offset for plotting multiple IRs
%% 

T60 = zeros(2, 2);
leave_out_first = round(0.05*fs);  %leave out first few ms before calculating envelope
smbeta = 50;  %smoothing window for energy envelope calculation (in ms)

for j = 1:2
            
    ir = h(leave_out_first:end, j); %leave out first 50 ms;
    ir_env = energy_envelope(ir,fs,smbeta);
    t_env = linspace(0,length(ir)/fs,length(ir_env)).';
    [taus, ~, ir_fit] = fit_two_stage_decay(ir_env, t_env, [rt60_pi, rt60_0]);
    T60(j, :) = taus * log(1000);
    c1 = log(ir_fit(1));     %calculate the right y-intercept assuming slope is right
    line1 = c1 - t_env./taus(1) ;
    c2 = log(ir_fit(end)) + t_env(end)/taus(2);
    line2 = c2 - t_env./taus(2) ;
    tp_x = (c2-c1)/(1/taus(2) - 1/taus(1));
    tp_y = c1 - tp_x/taus(1);


    %% plot 2 stage decay with and without FFM
    figure(j+101);
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
ylim([-120, 180]);xlim([0 L]); 
xlabel('Time (s)');
ylabel('Magnitude (dB)');

end 
