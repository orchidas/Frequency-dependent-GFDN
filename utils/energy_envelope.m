function [ env ] = energy_envelope( sig, fs, beta )
%beta is the amount of smoothing (in ms)

    len = length(sig);
    staps = round(beta/2 * fs/1000);
    bs = hanning(2*staps-1) / sum(hanning(2*staps-1));

    env = real(sqrt(fftfilt(bs, [zeros(staps,1); sig.^2; zeros(staps,1)])));
    env = env(2*staps-1+(1:len),:);


end
