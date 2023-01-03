function [P,T, signComb] = make_more_paraunitary(A)
%%
% AA is an NxN coupling matrix with coupling filters on the off-diagonals
% and zeros on the diagonals. P is a (somewhat) paraunitary version of A
% found by solving the Procrustes problem at each frequency bin

nRooms = size(A,1);

% frequency-response
AA = fft(ifftshift(A,3),[],3); %% is real-valued

% energy per column must be lower than 1, else, there is no way to create a
% orthogonal matrix. The same is true of row, but that's covered by
% symmetry.
E = squeeze(sum(AA.^2,1));

% create missing diagonal from remaining energy 
DD = zeros(size(AA));
for it = 1:size(AA,3)
    DD(:,:,it) = diag(sqrt(1 - E(:,it)));
end

squeeze(sum((AA+DD).^2,1)) % energy matches

% There a few sign options. Let's check which one is the best.
signCombinations = logical_indices_for_power_set(nRooms*(nRooms-1));

MM = AA+DD;
frobError = zeros(1,length(signCombinations));
for it = 1:length(signCombinations)
   signComb = signCombinations{it}*2 - 1;
   
   % set sign
   TT = signFlipVersion(MM,signComb);
   
   % solve Procrustes per-frequency
   PP = ProcrustesPerBin(TT);
   % Frobenius norm
   frobError(it) = sum(sum(sum(abs(PP - TT).^2)));
   
end

[~, bestSignCombo] = min(frobError);

signComb = signCombinations{bestSignCombo}*2 - 1;
% set sign
TT = signFlipVersion(MM,-signComb); 


% before solving Procrustes per bin
T = fftshift(ifft(TT,[],3),3);

% solve Procrustes per-frequency
PP = ProcrustesPerBin(TT);
   
% go back to time-domain
P = fftshift(ifft(PP,[],3),3);



%% helper functions
function TT = signFlipVersion(TT,signComb)
   % flip diagonal signs
   TT(1,1,:) = TT(1,1,:) * signComb(1);
   TT(2,2,:) = TT(2,2,:) * signComb(2);
   TT(3,3,:) = TT(3,3,:) * signComb(3);
   
   % flip off diagonal signs
   TT(2,1,:) = TT(2,1,:) * signComb(4);
   TT(1,2,:) = TT(1,2,:) * signComb(4);
   
   TT(2,3,:) = TT(2,3,:) * signComb(5);
   TT(3,2,:) = TT(3,2,:) * signComb(5);
   
   TT(3,1,:) = TT(3,1,:) * signComb(6);
   TT(1,3,:) = TT(1,3,:) * signComb(6);
end
  
function PP = svdPerBin(TT)

    PP = zeros(size(TT,[1 3]));
    for f = 1:size(TT,3)
       s = (svd(TT(:,:,f)));  
       PP(:,f) = s;
    end
end


function PP = ProcrustesPerBin(TT)

    PP = zeros(size(TT));
    for f = 1:size(TT,3)
       PP(:,:,f) = Procrustes(TT(:,:,f));
    end
end
   
function P = Procrustes(A)
    [U,S,V] = svd(A);
    P = U * V';
end

end