function [data_transmit NoCarriers] = Transmitter(data, NoPilots, trellis)
%%Transmitter of the OFDM system

%% Convolutionally encoding data
codedata = convenc(data, trellis);

%%Interleaving coded data
QAMbit = 4; %16-QAM: 2^4=16
matrix = reshape(codedata,size(codedata,2)/QAMbit, QAMbit); % codedata is a column vector
intlvddata = matintrlv(matrix',2,QAMbit/2)'; %Interleave
intlvddata = intlvddata';

%% Binary to decimal conversion
dec = bi2de(intlvddata','left-msb');

%% 16-QAM Modulation/16-PSK Modulation

y = qammod(dec,16);
%scatterplot(y);
%hModulator = comm.PSKModulator(M); % BitInput, true);
%hModulator.PhaseOffset = pi/16; %change the pashe offset to pi/16;
%y = step(hModulator, dec);
% constallation(y);

%% Pilot insertion
% pilt_data = PilotInsertion(y,NoPilots);

%% IFFT
Nocarriers = length(y);
ifft_sig = ifft(y',Nocarriers);

%% ADD Cyclic Extension
data_transmit = CyclicPrefixAdd(ifft_sig);