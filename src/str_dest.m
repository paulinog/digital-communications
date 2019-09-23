function text_output = str_dest(array_bin)
%STRING_DESTINATION Recebe um vetor binario e converte em texto

len_src = length(array_bin);

output_bin = reshape(array_bin, [len_src/8,8]);
text_output = sprintf('%s',bi2de(output_bin));

end
