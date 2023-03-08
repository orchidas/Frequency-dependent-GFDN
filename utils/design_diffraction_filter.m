function [b_diff, a_diff , win_len, tau_avg] = design_diffraction_filter(aperture, Fs, c, varargin)
%% Transmission coefficient through aperture connecting two rooms
% INPUTS : 
% aperture - aperture size, 
% fs - sampling rate, 
% c - speed of sound, 
% OPTIONAL :
% des_type - design type. exact or approximate
% OUTPUTS : 
% b_diff - FIR filter coefficients
% a_diff - 1
% win_len - length of the filter
% OUTPUTS : 
% b_diff, a_diff - diffraction filter coefficients
% win_len - length of the filter
% tau_avg - theoretical filter in the frequency domain

    if nargin > 3
        des_type = varargin{1};
    else
        des_type = 'approx';
        tau_avg = [];
    end
    
    order = 4*(aperture * Fs/c);
    sigma =  pi* aperture^2;     %transmission cross-section, assuming circular aperture
    win_len = round(order);
    
    % make sure window length is odd
    if (mod(win_len,2) == 0)
        win_len = win_len + 1;
    end
    
    %% design choice 1 - correct way
    if strcmp(des_type, 'exact') || strcmp(des_type, 'pm')
    
        nfreq = 2^12 ; % 2^nextpow2(win_len);
        freqs = linspace(0, Fs/2, nfreq);
        k = 2*pi*freqs/c;   %wave number
    
        % get average transmission by integrating over all angles of incidence
        fun = @(angle) sigma * (sinc((k * aperture)/pi * sin(angle)).^2).*sin(angle); 
        tau_avg = 2*integral(fun, 0 , pi/2, 'ArrayValued',true); 

        % we know that the transmission coefficient should be 1 for high frequencies
        cutoff = find(k*aperture >= 5.0,1,'first');
        tau_avg(cutoff:end - 1) = 1;
        tau_avg(end) = 0;
        startFrom = find(k*aperture >= 3.0, 1, 'first');
        %interpolation to avoid abrupt jump
        tau_avg(startFrom:cutoff) = linspace(tau_avg(startFrom), tau_avg(cutoff), length(startFrom:cutoff));

        if strcmp(des_type,'pm')
            %use parks mcClellan to design FIR filter    
            b_diff = firpm(win_len - 1, freqs/(Fs/2), tau_avg).'; 
            win_len = length(b_diff);
            a_diff = 1;
        else
            % inverse FFT from integral - this is what filter should be
            tau_spec = [tau_avg, conj(tau_avg(2:end-1))];
            % a real signal will have a hermitian symmetric spectrum
            tau_filt = fftshift(ifft(tau_spec, 'symmetric'));
        
            b_diff = tau_filt(nfreq-(win_len-1)/2 : nfreq+(win_len-1)/2).';  
            a_diff = 1;
        end
    
    %% design choice 2 - approximate sinc^2 response with bartlett window,
    % faster but less accurate
    else
           
        win = bartlett(win_len);
        a_diff = 1;
        %normalize the window to have proper magnitude
        b_diff = win.* (Fs*sigma/(4*c)*(pi/2)/win_len);     
    end


end