function codeData = symbolConverter(message, originalSymbolSet, targetSymbolSet)

numberOfElements = length(message);

codeDataTemp = [];

for i = 1:numberOfElements
    if message(i) == originalSymbolSet(1)
        codeDataTemp = [codeDataTemp targetSymbolSet(1)];
    else
        codeDataTemp = [codeDataTemp targetSymbolSet(2)];
    end
end

codeData = codeDataTemp;