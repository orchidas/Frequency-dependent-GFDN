function [X, ax] = ftgram(x, fs, typename, varargin)
% FBGRAM - compute, plot short-time filter bank output
%
% [X, AX] = fbgram(x, FS, TYPE, VARARGIN) returns X, the short-time filter bank
% output of the input x, computed using stfb().  The STFB is plotted using
% the sampling rate FS in Hz according to the string TYPE and parameters
% specified in name, value pairs.  The variable TYPE can be 'rir' or 'music';
% it sets the defaults for the spectrogram and plotting parameters described
% below.
%
% NAME = [RIR_DEFAULT MUSIC_DEFAULT];   %% DESCRIPTION, UNITS
%
% STFB parameters
% 'nbins' = [96 480];  %% filter bank window half length, bins
% 'nskip' = nbins/2;  %% hop size, samples
% fb = 125*2.^([-2:0.5:7]');  %% band center frequencies, Hz
% order = 2;  %% filter bank band splitting filter order, poles
% 
% spectrogram image axes
% 'dbrange' = [80 60];   %% gram dynamic range, dB
% 'logf' = [true true];  %% logarithmic frequency axis, indicator
% 'logt' = [true false]; %% logarithmic time axis, indicator
% 'ms' = [true false];   %% time/frequency axis ms/kHz scaling, indicator
% 
% waveform onset trimming
% 'trim' = [true false]; %% trim waveform onset, indicator
% 'preroll' = 10; %% onset zeropad duration, milliseconds
% 'onsetlevel' = 1e-2;   %% onset level, fraction
% 
% waveform plot parameters
% 'waveform' = [true false];    %% plot waveform, indicator
% 'tanhflag' = [true false];    %% hyperbolic tangent saturation, indicator
% 'tanhbeta' = 5; %% hyperbolic tangent saturation parameter, ratio
%
% The return variable AX contains the plot axes, one for each spectrogram
% channel, the first being the waveform plot, if present.
%
% See Also: STFT, IRGRAM.

% Created: 18-Mar-2012.
% Revised: 18-Mar-2012, JSA, v1.
% Version: v1.


%% initialization, parse input

% initialize type defaults
nbins_default = [96 480];    %% dft half length, bins
fb_default = 125*2.^([-2:0.5:7]');  %% band center frequencies, Hz
order_default = 2;  %% filter bank band splitting filter order, poles

dbrange_default = [80 60 60];   %% gram dynamic range, dB
logf_default = [true true false];   %% logarithmic frequency axis, indicator
logt_default = [true false false];  %% logarithmic time axis, indicator
logtmin_default = [10 10 10];  %% logarithmic time axis offset, milliseconds
ms_default = [true false false];    %% time/frequency axis ms/kHz scaling, indicator

trim_default = [true false false];  %% trim waveform onset, indicator
preroll_default = 10*[1 1 1];   %% onset zeropad duration, milliseconds
onsetlevel_default = 1e-2*[1 1 1];  %% onset level, fraction

waveform_default = [true false false];  %% plot waveform, indicator
tanhflag_default = [true false false];  %% hyperbolic tangent saturation, indicator
tanhbeta_default = 2*[1 1 1];   %% hyperbolic tangent saturation parameter, ratio

% set signal type
switch typename,
    case 'rir',
        % room impulse response
        type = 1;
    case 'music',
        % music input
        type = 2;
    otherwise,
        % music default
        type = 2;
end;

% parse input
p = inputParser;

% waveform, sampling rate
p.addRequired('x', @(x)isnumeric(x));   %% waveform, signal matrix
p.addRequired('fs', @(x)isnumeric(x) && x>0);   %% sampling rate, Hz
p.addRequired('typename', @(x)ischar(x));   %% sampling rate, Hz

% stft parameters
p.addParamValue('nbins', nbins_default(type), @(x)isnumeric(x));
p.addParamValue('nskip', 0, @(x)isnumeric(x));
p.addParamValue('fb', fb_default, @(x)isnumeric(x));
p.addParamValue('order', order_default, @(x)isnumeric(x));

% spectrogram image axes
p.addParamValue('dbrange', dbrange_default(type), @(x)isnumeric(x));
p.addParamValue('logf', logf_default(type), @(x)islogical(x));
p.addParamValue('logt', logt_default(type), @(x)islogical(x));
p.addParamValue('logtmin', logtmin_default(type), @(x)isnumeric(x));
p.addParamValue('ms', ms_default(type), @(x)islogical(x));

% waveform onset trimming
p.addParamValue('trim', trim_default(type), @(x)islogical(x));
p.addParamValue('preroll', preroll_default(type), @(x)isnumeric(x));
p.addParamValue('onsetlevel', onsetlevel_default(type), @(x)isnumeric(x));

% waveform plot parameters
p.addParamValue('waveform', waveform_default(type), @(x)islogical(x));
p.addParamValue('tanhflag', tanhflag_default(type), @(x)islogical(x));
p.addParamValue('tanhbeta', tanhbeta_default(type), @(x)isnumeric(x));

p.parse(x, fs, typename, varargin{:});

% assign variables
nbins = p.Results.nbins;
if (p.Results.nskip > 0);
    nskip = p.Results.nskip;
else,
    nskip = nbins/2;
end;
fb = p.Results.fb(find(p.Results.fb < fs/2));
order = p.Results.order;

dbrange = p.Results.dbrange;
logf = p.Results.logf;
logt = p.Results.logt;
logtmin = p.Results.logtmin;
ms = p.Results.ms;

trim = p.Results.trim;
preroll = p.Results.preroll;
onsetlevel = p.Results.onsetlevel;

waveform = p.Results.waveform;
tanhflag = p.Results.tanhflag;
beta = p.Results.tanhbeta;


%% condition input, form spectrogram

% find input signal size
nsamp = size(x,1);
if (nsamp == 1),
    % make x a column
    x = x(:);
    [nsamp, nc] = size(x);
end;
[nsamp, nc] = size(x);  %% signal length, samples; channel count, channels

% trim waveform onset
if trim,
    istart = find(mean(abs(x)/max(max(abs(x))),2) > onsetlevel, 1) - round(preroll*fs/1000);
    x = x(max(1,istart):end,:);
end;
nsamp = size(x,1);

% compute, normalize spectrogram
X = stfb(x, nbins, nskip, fb*2/fs, order);
nframes = size(X,2)/nc;

y = x/max(max(abs(x)));
Y = 20*log10(abs(X)/max(max(abs(X)))+eps);


%% plot waveform

if waveform,
    % define axis
    figure(gcf);
    ax = subplot(max(nc+1,3),1,1);

    % tanh scaling
    if tanhflag,
        y = tanh(beta*y)/tanh(beta);
    end;

    % define time axis
    t = [0:nsamp-1]/fs;
    tscale = (~ms)*1 + ms*1000;

    % plot waveform, label time axes
    if logt,
        semilogx(tscale*(t+logtmin/1000), y); grid;
        xlim(tscale*[logtmin/1000 logtmin/1000+nsamp/fs]);
        n = get(gca,'Xtick');
        set(gca,'XTickLabel',sprintf('%g |', n'));
    else,
        plot(tscale*t, y); grid;
        xlim(tscale*[0 nsamp/fs]);
    end;

    % y-axis labels
    if tanhflag,
        ylabel('Amplitude (dB)');
        set(gca,'YTick', sort(kron([-1 1], tanh(10.^-([0:10:30]/20)*beta))));
        set(gca,'YTickLabel',['  0'; '-10'; '-20'; '-30'; '-30'; '-20'; '-10'; '  0'])
        ylim(tanh(10^(1/20)*beta)*[-1 1]);
    else,
        ylim([-1 1]);
        ylabel('Amplitude');
    end;

end;


%% plot spectrogram

% define time, frequency axes
tscale = (~ms)*1 + ms*1000;

np = ceil(nbins/(2*nskip));
nq = ceil((nsamp-(nframes-1)*nskip-nbins/2)/nskip);
t = tscale * ([-np:nframes-1+nq]*nskip + nbins/2)/fs;
f = 1/tscale * [eps; exp(mean(log([fb(1:end-1) fb(2:end)]), 2)); fs/2];
Y = Y([[1:end] end],:);

% loop through specgtrograms
for s = [1:nc],

    % get axes
    if (nc == 1),
        if waveform,
            ax(2) = subplot(3,1,[2 3]);
        else,
            ax(1) = subplot(1,1,1);
        end;
    else,
        ax(waveform+s) = subplot(nc+waveform,1,s+waveform);
    end;

    % display spectrogram
    offset = logt * tscale*logtmin/1000;
    surf(t+offset, f, Y(:, (s-1)*nframes + [ones(1,np) [1:nframes] nframes*ones(1,nq)]), 'edgecolor', 'none');
    axis tight;
    view(0,90);

    % time, frequency scaling
    if ms,
        if (s == nc),
            xlabel('Time (ms)');
        end;
        ylabel('Frequency (kHz)');
    else,
        if (s == nc),
            xlabel('Time (s)');
        end;
        ylabel('Frequency (Hz)');
    end;

    % scale, label frequency axis
    if logf,
        set(gca, 'yScale', 'log');
    end;

    % scale, label time ax1s
    if logt,
        set(gca, 'xScale', 'log');

        xlim(tscale*[logtmin/1000 logtmin/1000+nsamp/fs]);
        n = get(gca,'Xtick');
        set(gca,'XTickLabel',sprintf('%g |', n'));

    end;

    % scale, label frequency axis
    if logf,
        divs = [50 100 200 500 1000 2000 5000 10000 20000]/tscale;
        set(gca, 'ytickmode', 'manual');
        set(gca, 'ytick', divs);

        ylim([min(divs) max(divs)]);

    else,
        ylim([0 20000]/tscale);

    end;

    % display color bar
    caxis([-dbrange 0])
    B = colorbar();
    if waveform*(nc == 1),
        set(B, 'Position', [0.916 0.11 0.015 0.517]);
    elseif (nc == 1),
        set(B, 'Position', [0.916 0.11 0.015 0.815]);
    else,
        temp = get(B, 'Position');
        set(B, 'Position', [0.916 temp(2)+0.004 0.015 temp(4)]);
    end;
    colormap(jet);
    ylabel(B,'Energy (dB)');

end;

% link x-axes of plots:
if (waveform + nc-1),
    linkaxes(ax,'x');
end;


end
