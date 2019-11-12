function a = mapper(array_bin, symbols_set)
%MAPPER Mapeador de bits para simbolos

%% OOK
% M=2;
% symbols_set = [0, 1];
%a = array_bin
%% M-PAM
% M=4;
% symbols_set = [-3, -1, 1, 3];

%% symbol mapper
M = length(symbols_set);
b = log2(M);
len_src = length(array_bin);

a_mapped = reshape(array_bin, [b, len_src/b])';
symbol_index = bi2de(a_mapped);
a = symbols_set(symbol_index+1);
end 
