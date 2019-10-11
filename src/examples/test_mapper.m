%% Exemplo - Projeto Final IE533
% gpaulino
clc; clearvars; close all;
disp('Final Project (IE533)')
%% Source
disp('Enter text:')
input_str = input('','s');
input_bin = de2bi(real(input_str),8);
array_bin = reshape(input_bin, [1,size(input_bin,1)*size(input_bin, 2)]);
len_src = length(array_bin);
%% mapper
% OOK
% M=2;
% symbols_set = [0, 1];
%a = array_bin

% M-PAM
M=4;
symbols_set = [-3, -1, 1, 3];

% symbol mapper
b=log2(M);
a_mapped = reshape(array_bin, [b, len_src/b])';
symbol_index = bi2de(a_mapped);
a = symbols_set(symbol_index +1);

%% modulate
% e.g. 4-PAM
% pam4
z = a;
% use sound-function

%% demapper
out_index = zeros(1,length(z));
for i = 1 : length(z)
    out_index(i) = find(symbols_set == z(i));
end
out_bin = de2bi(out_index-1);
out_array = reshape(out_bin', [1,len_src]);

%% Destination
output_bin = reshape(out_array, [len_src/8,8]);
text_output = sprintf('%s',bi2de(output_bin));

disp('-----')
disp('Output text:');
disp(text_output)