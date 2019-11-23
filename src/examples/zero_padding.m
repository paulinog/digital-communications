function [num_padding] = zero_padding(num_bits, N)
%ZERO_PADDING
num_multiple = ceil(num_bits/N)*N;
num_padding = num_multiple - num_bits;
end

