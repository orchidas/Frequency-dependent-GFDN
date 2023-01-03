function [zAbsorption, feedbackMatrix, inputGain, outputGain, directGain] = convert_to_FDNTB_params(nGrp,nSize,b,a,theta,b_drive,c_drive,d,varargin)
%% Convert FDN parameters to FDNTB usable forms
% INPUTS:
% nGrp - number of groups in GFDN
% nSize - size of each group (vector)
% b,a - IIR absorption filter coefficients
% theta - mixing angles (for mixing matrix)
% b_drive, c - input and output gains
% d - direct path gain
% OPTIONAL:
% couplingMatrix - coupling matrix among groups
% OUTPUTS:
% zAbsorption - absorption filters as diagonal matrix of filters
% feedbackMatrix - mixing matrix
% inputGain, outputGain - input and output gain matrices
% directGain - direct gain vector

p = inputParser;
p.KeepUnmatched = true;
addParameter(p,'couplingMatrix',[]);
parse(p,varargin{:});
couplingMatrix = p.Results.couplingMatrix;


% input and output coefficients
inputGain = b_drive;
outputGain = blkdiag(c_drive.');
directGain = d*ones(nGrp,1);


% T60 filters
nDel = sum(nSize);
absorption.b = reshape(b, [nDel,1,size(b,2)]);
absorption.a = reshape(a, [nDel,1,size(a,2)]);
zAbsorption = zTF(absorption.b,absorption.a,'isDiagonal',true);



if nGrp > 1
    % block mixing matrix
    M_block = get_block_mixing_matrix(nSize, nGrp, theta);
    % feedback matrix
    degree = size(couplingMatrix,3);
    feedbackMatrix = zeros(nDel, nDel, degree);
    for k = 1:degree
        feedbackMatrix(:,:,k) =  kron(couplingMatrix(:,:,k), ones(nDel/nGrp, nDel/nGrp)) .* M_block;
    end
else
   feedbackMatrix = get_mixing_matrix(nSize, theta);
end





end