function [A] = two_room_coupling_matrix(aperture, b_diff, areas, varargin)

%% Coupling matrix for two coupled rooms
% INPUTS:
% aperture : aperture radius (in metre)
% b_diff : FIR diiffraction filter coefficients (vector)
% areas : areas of 2 rooms (vector)
% OPTIONAL:
% asymmetry : 0 or 1, whether to have asymmetric coupling (boolean)
% OUTPUTS:
% A : NxNxL FIR polynomial matrix

p = inputParser;
p.KeepUnmatched = true;
addParameter(p, 'asymmetry', 1);
parse(p,varargin{:});
asymmetry = p.Results.asymmetry;
addpath('fdnToolbox/.');

S1 = areas(1); S2 = areas(2);
N = 2;

% Create some coupling filters with Bartlett windows
couplingLength = length(b_diff); % should be all odd
% add coupling gains
couplingGain = (pi*aperture^2)/sqrt(S1*S2);

% build up coupling filter matrix
A = zeros(N,N,couplingLength);
for i = 1:N
    for j = i+1:N
        w = b_diff;
        g = couplingGain;
        z = zeros(couplingLength,1);
        shift = floor(length(w)/2);
        z((-shift:shift)+ceil(couplingLength/2)) = w;
        A(i,j,:) = g*flip(z);
        A(j,i,:) = -g*z; % for symmetry
    end
end


alpha = squeeze(A(1,2,:));
aa = conv(alpha, flip(alpha));
aa1 = aa;
% aa1(ceil(end/2)) = 1-aa1(ceil(end/2));
% aa1 = -aa1;

% find roots and invert those within unit circle
% ra = roots(aa1);
% rx = ra(abs(ra) < 1);
% x = poly(rx).';
% xx = conv(x,flip(x));
% f = xx(3) / aa(3);
% x = x/sqrt(abs(f));

% just a linear phase FIR filter
X_mag = sqrt(1 - abs(fft(aa1)));
F = linspace(0,1,couplingLength+1);
x = firls(couplingLength-1, F, X_mag(1:couplingLength+1));


A(1,1,:) = x;
A(2,2,:) = flip(x);

% [isP, ~, ~] = isParaunitary(A);
% assert (isP == 1, 'Coupling matrix not paraunitary');

% multiply with scaling factor for asymmetric coupling
if asymmetry
    d = sqrt(S2/S1);
    A(1,2,:) = d * A(1,2,:);
    A(2,1,:) = 1/d * A(2,1,:);
end



end

