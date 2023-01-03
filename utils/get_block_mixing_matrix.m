function [M_block] = get_block_mixing_matrix(nSize, nGrp, theta)

%% Returns an orthonormal block mixing matrix (scalar) for a GFDN
% INPUTS:
% nSize - number of delay lines in each group
% nGrp - number of groups
% theta - the vector of mixing angles
% OUTPUTS:
% M_block - block mixing matrix (nSizexnGrp x nSizexnGrp)
%%

M_block = zeros(sum(nSize),sum(nSize));
row = 1; 

for i = 1:nGrp
    col = 1;
    for j = 1:nGrp
        if(i == j)
            M_block(row : row+nSize(i)-1, col:col+nSize(j)-1) = get_mixing_matrix(nSize(i), theta(i));
        else
            M_block(row : row+nSize(i)-1, col:col+nSize(j)-1) = (get_mixing_matrix(nSize(i), theta(i)/2) * ...
                get_mixing_matrix(nSize(j), theta(j)/2));
        end
        col = col + nSize(j);
    end
    row = row + nSize(i);
end

