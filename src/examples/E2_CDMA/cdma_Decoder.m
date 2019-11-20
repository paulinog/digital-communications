function recoveredMessage = cdma_Decoder(receivedMessage, spreadingCode, symbolsSet, transmissionSymbolSet)

codeLen = length(spreadingCode);
receivedMessageDataLen = length(receivedMessage)/codeLen;

multiplicationSpreading = [];
for i = 1:receivedMessageDataLen
    multiplicationSpreading = [multiplicationSpreading symbolConverter(spreadingCode, symbolsSet, transmissionSymbolSet)];
end

multiplicationSpreading;
multipliedMessage = multiplicationSpreading .* receivedMessage;

recoveryMessageTemp = [];
for i = 1:receivedMessageDataLen
    recoveryMessageTemp = [recoveryMessageTemp sum(multipliedMessage((1:codeLen) + (i-1)*codeLen))];
end

% recoveryMessageTemp = recoveryMessageTemp/codeLen
% recoveryMessageTemp = symbolConverter(recoveryMessageTemp, transmissionSymbolSet, symbolsSet);

%% SLICER
% recoveryMessageTemp = recoveryMessageTemp/codeLen
recoveryMessageTemp;

for i = 1:length(recoveryMessageTemp)
    if recoveryMessageTemp(i) <= 0
        recoveryMessageTemp(i) = 0;
    else
        recoveryMessageTemp(i) = 1;
    end
end

recoveredMessage = recoveryMessageTemp;
% recoveredMessage = [1 1]