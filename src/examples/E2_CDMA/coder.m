function codeData = coder(data, spreadingCode)

numBits = length(data);
codeDataTemp = [];

for i = 1:numBits
    %codeDataTemp = [codeDataTemp xor(data(i), spreadingCode)];
    codeDataTemp = [codeDataTemp data(i)*spreadingCode];
end

codeData = codeDataTemp;