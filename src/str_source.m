function [array_bin, len_src] = str_source(text_input)
%STRING_SOURCE Recebe um texto digitado pelo usuario e converte em vetor
% binario

input_bin = de2bi(real(text_input),8);

array_bin = reshape(input_bin, [1,size(input_bin,1)*size(input_bin, 2)]);
len_src = length(array_bin);

end
