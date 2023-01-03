function[M] = get_mixing_matrix(N, theta, type)
%% Get unitary mixing matrix of FDN with N delay lines by constructing 
% rotation or reflection matrices with angle theta
% INPUTS :
% N - order of matrix (int)
% theta - mixing angle (radians)
% type - rotation or reflection matrix
% OUTPUTS :
% M - NxN mixing matrix 

    if nargin < 3
        type = 'rot';
    end
    
    %rotation or reflection matrices
    rot_2 = [cos(theta) -sin(theta);sin(theta) cos(theta)];
    ref_2 = [cos(2*theta) sin(2*theta);sin(2*theta) -cos(2*theta)];
    if strcmp(type,'rot')
        M = rot_2; 
    else
        M = ref_2;
    end
    
    n = log2(N);
    for k = 1:n-1
        M = kron(rot_2,M);
    end
end