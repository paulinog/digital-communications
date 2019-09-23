function out_array = demapper(z, symbols_set)
%DEMAPPER Demapeador de simbolos para bits

M = length(symbols_set);
b=log2(M);
len_src = length(z);
out_index = zeros(1, len_src);

for i = 1 : len_src
    out_index(i) = find(symbols_set == z(i));
end

out_bin = de2bi(out_index-1, b);
out_array = reshape(out_bin', [1, numel(out_bin)]);

end 
