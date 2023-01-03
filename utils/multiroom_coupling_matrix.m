function [couplingMatrix, S] = multiroom_coupling_matrix(Nrooms, filter_coeffs, aperture, couplingAreas)

%% Polynomial paraunitary coupling matrix for N > 2 coupled rooms
% INPUTS:
% Nrooms - number of coupled rooms
% filter_coeffs - Nrooms x Nrooms x degree diffraction filter coefficient
% matrix
% aperture - matrix containing relative aperture sizes
% couplingAreas - matrix containing the coupling areas for each room pair
% OUTPUTS:
% couplingMatrix - Nrooms x Nrooms x degree coupling matrix
% S - per bin Procrustes solution of couplingMatrix


win_len = zeros(Nrooms, Nrooms);
for i = 1:Nrooms
    for j = 1:Nrooms
        win_len(i,j) = find(filter_coeffs(i,j,:) == 0,1);
        asymmetric_coupling_coeff = sqrt(couplingAreas(i,j)/couplingAreas(j,i));
    end
end

degree = max(max(win_len))-1;
couplingMatrix = zeros(Nrooms, Nrooms, degree);
zeroPos = (degree+1)/2;     %position of 0th order coefficient


for i = 1:Nrooms
    for j = 1:Nrooms
        if ( i == j)
            continue;
        else
            shift = (win_len(i,j)+1)/2;
            couplingMatrix(i,j,zeroPos-shift+1:zeroPos+shift-1) = (pi*(aperture(i,j)^2))/sqrt(couplingAreas(i,j)*couplingAreas(j,i))...
                *filter_coeffs(i,j,1:win_len(i,j));
        end
    end
end


% make coupling matrix more paraunitary
couplingMatrix = make_more_paraunitary(couplingMatrix);
S = svdPerBin(fft(couplingMatrix,2*degree,3));


% add asymmetric coupling
E = diag(asymmetric_coupling_coeff);
Einv = diag(1./asymmetric_coupling_coeff);
couplingMatrix_scaled = zeros(size(couplingMatrix));

for k = 1:degree+1
    couplingMatrix_scaled(:,:,k) = E * couplingMatrix(:,:,k) * Einv;
end



function PP = svdPerBin(TT)

    PP = zeros(size(TT,[1 3]));
    for f = 1:size(TT,3)
       s = (svd(TT(:,:,f)));  
       PP(:,f) = s;
    end
end


end