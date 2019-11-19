%% Codificacao e decodificacao CDMA
clear all;
clc;

%% Definicao dos dados da mensagem
symbolsSet = [0 1];
transmissionSymbolSet = [-1 1];
numberOfUsers = 1;

data = [];

% Data for tests
% data = [[0 0];[1 0];[1 1]];

% Data generated randomly
% dataLen = 200;
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
codeData = [];

for user = 1:numberOfUsers
    codeData = [codeData; coder(dataConverted(user,:), spreadingCodeConverted(user,:))];
end

codeData;

%% Composed Message

if numberOfUsers == 1
    transmitedMessage = codeData;
else
    transmitedMessage = sum(codeData);
end

transmitedMessage;

%% Modulation

%% Transmission - Tx
s = transmitedMessage;

%% Channel
r = s; %Bypass

%% Receiver - Rx
receivedMessage = r;

%% Decoder
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
recoveredMessage ~= data;
err = sum(sum(recoveredMessage ~= data));
if (err > 0)
    error([num2str(err) ' errors'])
else
    disp('no errors')
end
