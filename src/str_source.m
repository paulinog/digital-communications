function [array_bin, len_src] = str_source()
%STRING SOURCE Receives a text typed by the user and converts it to a
% binary array

input_str = input('','s');
input_bin = de2bi(real(input_str),8);

array_bin = reshape(input_bin, [1,size(input_bin,1)*size(input_bin, 2)]);
len_src = length(array_bin);

end
