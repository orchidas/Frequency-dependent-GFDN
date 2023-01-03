function [] = ftgram2N3S( x, fs, figtitle, saveflag )
% yields plot with nice labels
% withthree, 1/3 spectrograms

tempfig = figure;
[~, ax] = ftgram(x, fs, 'rir',  'waveform', false, 'dbrange', 60, 'tanhflag', false, 'ms', false, 'trim', false, 'logt', false, 'nbins', 1024, 'nskip', 16);

fig = figure('Units','inches', 'Position',[0 0 3.29 2.5],'PaperPositionMode','auto');
    set(gca, 'FontUnits','points', 'FontWeight','normal', 'FontSize',8, 'FontName','Times');

nax(1) = subplot(3,1,1);
    set(gca, 'FontUnits','points', 'FontWeight','normal', 'FontSize',8, 'FontName','Times');

nax(2) = subplot(3,1,2);
    set(gca, 'FontUnits','points', 'FontWeight','normal', 'FontSize',8, 'FontName','Times');

nax(3) = subplot(3,1,3);
    set(gca, 'FontUnits','points', 'FontWeight','normal', 'FontSize',8, 'FontName','Times');


sg1 = get(ax(1), 'Children');
copyobj(sg1, nax(1));
sg2 = get(ax(2), 'Children');
copyobj(sg2, nax(2));
sg3 = get(ax(3), 'Children');
copyobj(sg3, nax(3));

close(tempfig);

linkaxes(nax, 'x');


axes(nax(1));
axis tight;
set(gca, 'YScale', 'log');
ylim([20 20000]);
% xlabel('Time (s)');
% yl2 = ylabel('Frequency (kHz)');
% yl2p = get(yl2, 'position');
yticks([100, 1000, 10000]);
yticklabels({'0.1', '1', '10'});
xticklabels({});
cb = colorbar;
caxis([-60 0]);
nax1p = get(nax(1), 'Position');
nax1p(3) = 0.7; % adjust width to compenstate for color bar
nax1p(2) = 0.7; % move upwards to make room for x label
nax1p(4) = 0.2; % make subplots slightly larger

cbl = get(cb,'Position');
cbl(1) = 0.84;
cbl(2) = 0.2;
cbl(3) = 0.03;
set(cb, 'Position', cbl);
set(nax(1), 'Position', nax1p);
ylabel(cb,'Energy (dB)');
    set(gca, 'FontUnits','points', 'FontWeight','normal', 'FontSize',8, 'FontName','Times');
    
    
axes(nax(2));
axis tight;
set(gca, 'YScale', 'log');
ylim([20 20000]);
% xlabel('Time (s)');
% yl3 = ylabel('Frequency (kHz)');
% yl3p = get(yl3, 'position');
yticks([100, 1000, 10000]);
yticklabels({'0.1', '1', '10'});
xticklabels({});

yl2 = ylabel('Frequency (kHz)');
yl2p = get(yl2, 'position');

% cb2 = colorbar;
caxis([-60 0]);
nax2p = get(nax(2), 'Position');
nax2p(3) = 0.7; % adjust width to compenstate for color bar
nax2p(2) = 0.45; % move upwards to make room for x label
nax2p(4) = 0.2; % make subplots slightly larger

set(nax(2), 'Position', nax2p);


axes(nax(3));
axis tight;
set(gca, 'YScale', 'log');
ylim([20 20000]);
% xlabel('Time (s)');
% yl3 = ylabel('Frequency (kHz)');
% yl3p = get(yl3, 'position');
yticks([100, 1000, 10000]);
yticklabels({'0.1', '1', '10'});
% cb2 = colorbar;
caxis([-60 0]);
nax3p = get(nax(3), 'Position');
nax3p(3) = 0.7; % adjust width to compenstate for color bar
nax3p(2) = 0.2; % move upwards to make room for x label
nax3p(4) = 0.2; % make subplots slightly larger

xlabel('Time (s)');
    set(gca, 'FontUnits','points', 'FontWeight','normal', 'FontSize',8, 'FontName','Times');

% cbl2 = get(cb2,'Position');
% cbl2(1) = 0.84;
% cbl2(2) = 0.2;
% cbl2(3) = 0.03;
% set(cb2, 'Position', cbl2);
set(nax(3), 'Position', nax3p);
% ylabel(cb2,'Energy (dB)');
%     set(gca, 'FontUnits','points', 'FontWeight','normal', 'FontSize',8, 'FontName','Times');

% get positions of top and bottom spectrograms ylabel and colorbar
% yl2p = get(yl2, 'position');
% set(yl2,'position', [yl2p(1) 20 yl2p(3)]);
cbl = get(cb,'Position');
set(cb, 'Position', [cbl(1) cbl(2) cbl(3) cbl(4)*4.45]);


% axes(nax(1));
% set(yl1, 'position', [yl2p(1) yl1p(2:3)]);

checkIfTitleDefined = exist('figtitle', 'var');
if checkIfTitleDefined ~= 1
    figtitle = timecode;
end

drawnow;


if saveflag
audiowrite([figtitle(~isspace(figtitle)) '.wav'], 0.9*x/max(abs(x)), fs);
print([figtitle(~isspace(figtitle)) '.png'], '-dpng');
print([figtitle(~isspace(figtitle)) '.eps'], '-depsc');
end


end
