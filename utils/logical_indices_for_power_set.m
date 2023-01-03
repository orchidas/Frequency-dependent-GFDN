function logical_indices = logical_indices_for_power_set(set_size)

%% LOGICAL_INDICES_FOR_POWER_SET creates a cell array containing all possible
%combinations of logical indices to extract the members of a power set
%by indexing into the original set. This saves the memory overhead of
%storage of the set elements themselves, and instead only needs storage for
%the indices of the original set. 
%2^N cells of logical indices each of [1 x set_size].
%Converts the sequence of numbers from 0 to 2^set_size-1 to binary string
%representation. The characters of the string are turned into ascii
%characters. 48 is the ascii code for zero, so testing for equality will a 
%a matrix of logical indices, where the first row identifies all elements
%of the set.
logical_indices = mat2cell(double(dec2bin(0:2^set_size-1,set_size))==48,ones(2^set_size,1));