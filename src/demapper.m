function out_array = demapper(z, symbols_set)
%DEMAPPER

len_src = length(z);
out_index = zeros(1, len_src);

for i = 1 : len_src
    out_index(i) = find(symbols_set == z(i));
end

out_bin = de2bi(out_index-1);
out_array = reshape(out_bin', [1, numel(out_bin)]);

end 
