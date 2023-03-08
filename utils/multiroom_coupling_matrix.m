function [couplingMatrix, originalMatrix, degree] = multiroom_coupling_matrix(fs, c, Nrooms, aperture, absorp_coef, area, filter_type)

%% Polynomial paraunitary coupling matrix for N > 2 coupled rooms
% INPUTS:
% fs - sampling frequency (Hz)
% c - speed of sound in air (m/s)
% Nrooms - number of coupled rooms
% aperture - matrix containing relative aperture sizes
% absorp_coef - absorption coefficient of different rooms (Nrooms x 1)
% area - surface area of different rooms (Nrooms x 1)
% filter_type - diffraction filter type, exact or PM
% OUTPUTS:
% couplingMatrix - Nrooms x Nrooms x degree coupling matrix
% originalMatrix - Nrooms x Nrooms x degree original matrix (pre
% paraunitarization)
% degree - Degree of the polynomial coupling matrix

couplingAreas = ones(Nrooms, Nrooms);
filter_coeffs = zeros(Nrooms,Nrooms,1000);
win_len = zeros(Nrooms, Nrooms);

% get the diffraction filters
for i = 1:Nrooms
    for j = 1:Nrooms
        if i == j
            continue;
        else
            [b_diff, ~, win_len(i,j)] = design_diffraction_filter(aperture(i,j),fs,c,filter_type);
            filter_coeffs(i,j,1:win_len(i,j)) =  b_diff;
            
            % needed for asymmetric coupling coefficient
            couplingAreas(i,j) = absorp_coef(i).*area(i) + ...
                (pi*aperture(i,j)^2).*(1-absorp_coef(i));
        end
    end
end

% design the coupling matrix filters
degree = max(max(win_len));
couplingMatrix = zeros(Nrooms, Nrooms, degree);
zeroPos = (degree+1)/2;     %position of 0th order coefficient

for i = 1:Nrooms
    for j = 1:Nrooms
        if ( i == j)
            continue;
        else
            shift = (win_len(i,j)+1)/2;
            couplingMatrix(i,j,zeroPos-shift+1:zeroPos+shift-1) = ...
                (pi*(aperture(i,j)^2))/sqrt(couplingAreas(i,j)*couplingAreas(j,i))...
                *filter_coeffs(i,j,1:win_len(i,j));
        end
    end
end


% make coupling matrix more paraunitary
[couplingMatrix, originalMatrix, ~] = make_more_paraunitary(couplingMatrix);

% add asymmetric coupling
for i = 1:Nrooms
    for j = 1:Nrooms
        asymmetric_coupling_coeff = sqrt(couplingAreas(i,j)/couplingAreas(j,i));
    end
end

E = diag(asymmetric_coupling_coeff);
Einv = diag(1./asymmetric_coupling_coeff);
couplingMatrix_scaled = zeros(size(couplingMatrix));

for k = 1:degree
    couplingMatrix_scaled(:,:,k) = E * couplingMatrix(:,:,k) * Einv;
end



end