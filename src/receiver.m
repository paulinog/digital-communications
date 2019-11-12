%% Run Receiver
% Disciplina IE533A, 2s2019 (FEEC/UNICAMP)
% Guilherme Paulino, RA 117119
clc;
clearvars;
% close all;
disp('Final Project (IE533)')
disp('Binary FSK - RECEIVER')

%% Propriedades de Modulacao
symbols_set = [0 1];

fc0 = 440;
fc1 = 4*fc0;
fs = 4*fc1;

%% Codigo corretor de erros
trellis = poly2trellis(3, [5 7]);

%% Sincronismo de quadros (MLS)
k = 8;

%% Defina o numero de bits
numSymbol = input('Enter the number of bits to be received:');

%% Recebendo
numBits = 8;
numChannels = 1;
% audioRecordTime = 10; % em segundos
audio_rx = audiorecorder(fs, numBits, numChannels);

disp('Start recording audio')
record(audio_rx);
% pause(audioRecordTime);
input('Press ENTER to continue')

disp('Processing...')
stop(audio_rx)

r = getaudiodata(audio_rx);

%% Inicia um stopwatch timer
tic;

%% Demodulador
z = demod_fsk(r, fc0, trellis, k, numSymbol);

%% Demapeador
% output_bin = demapper(z, symbols_set);
output_bin = z;

%% Destino
output_str = str_dest(output_bin);
% Mostra a saida de texto
disp('-----')
disp('Output text:');
fprintf('%s', output_str)

%% Para o stopwatch timer
disp('-----')
disp('Total time to receive:');
toc;