function text_output = str_dest(array_bin)
%STRING DESTINATION Receive a binary array and converts it to a text output
% as concateneted strings

len_src = length(array_bin);

output_bin = reshape(array_bin, [len_src/8,8]);
text_output = sprintf('%s',bi2de(output_bin));

end
