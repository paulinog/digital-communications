%% Codificacao e decodificacao CDMA
clearvars;
clc;

%% Message Definition
symbolsSet = [0 1];
transmissionSymbolSet = [-1 1];
numberOfUsers = 1;

data = [];

% Data for tests
% data = [[0 0];[1 0];[1 1]];

% Data generated randomly
% dataLen = 2;
% for i = 1:numberOfUsers
%     data = [data; randsrc(1, dataLen, symbolsSet)]; % mensagem aleatoria
% end

% Data from a input
input_str = 'abcdefghijklmnopqrstuvxzwyABCDEFGHIJKLMNOPQRSTUVXZWY1234567890';
input_bin = str_source(input_str);
data = mapper(input_bin, symbolsSet);

% Binary conversion from high and low to (+/-)1
dataConverted = [];
for i = 1:numberOfUsers
    dataConverted = [dataConverted; symbolConverter(data(i,:), symbolsSet, transmissionSymbolSet)];
end

%% Spreading Code
codeLen = 10;
spreadingCode = [];

% Spreading Code for test
% spreadingCode = [[0 1 0 1]; [0 0 1 1]; [0 0 0 0]];

% Spreading Code generated randomly
for i = 1:numberOfUsers
    spreadingCode = [spreadingCode; randsrc(1, codeLen, symbolsSet)];
end

% Binary conversion from high and low to (+/-)1
spreadingCodeConverted = [];
for i = 1:numberOfUsers
    spreadingCodeConverted = [spreadingCodeConverted; symbolConverter(spreadingCode(i,:), symbolsSet, transmissionSymbolSet)];
end

%% VERIFICATION POINT
data;
spreadingCode;
dataConverted;
spreadingCodeConverted;

%% Spread Message
a = [];
for user = 1:numberOfUsers
    a = [a; coder(dataConverted(user,:), spreadingCodeConverted(user,:))];
end

a = symbolConverter(a, transmissionSymbolSet, symbolsSet);

%% Message encode
trellis = poly2trellis(3, [5 7]);
parity_ratio = 2;

k = 10;
sync_bits = 2^k-1;

a_enc = convenc(a, trellis);
syncVec = double(mls(k, 1) > 0);
a_sync = [syncVec a_enc syncVec]; % concatenate
% aSyncConverted = symbolConverter(a_sync, symbolsSet, transmissionSymbolSet);

%% Composed Message
if numberOfUsers == 1
    transmitedMessage = a_sync;
else
    transmitedMessage = sum(a_sync);
end

transmitedMessage;

%% Transmission - Tx
s = transmitedMessage;

figure()
stem(s)

%% Channel
r = s; %Bypass

%% Receiver - Rx
y_sync = r;

%% Self Correlation
self_corr = xcorr(y_sync, syncVec);
figure()
plot(self_corr);

max_peak_pos = (sync_bits/2)*0.95;

startFrame = sync_bits + min(find(self_corr(length(y_sync): end) > max_peak_pos));
endFrame = max(find(self_corr(length(y_sync): end) > max_peak_pos)) - 1;
frameLen = endFrame - startFrame;

%% Decoder
y_enc = y_sync(startFrame : endFrame);
y = vitdec(y_enc, trellis, 20,  'trunc', 'hard');

receivedMessage = y;

%% CDMA Decoder
recoveredMessage = [];
for i = 1:numberOfUsers
    recoveredMessage = [recoveredMessage; cdma_Decoder(receivedMessage, spreadingCode(i,:), symbolsSet, transmissionSymbolSet)];
end

%% Demapper
output_bin = demapper(recoveredMessage, symbolsSet);
output_str = str_dest(output_bin);

disp('Output text:');
fprintf('%s \n', output_str)

%% BER
err = sum(sum(recoveredMessage ~= data));
if (err > 0)
    error([num2str(err) ' errors'])
else
    disp('no errors')
end
