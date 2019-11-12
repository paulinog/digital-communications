%% Run Simulation
% Disciplina IE533A, 2s2019 (FEEC/UNICAMP)
% Guilherme Paulino, RA 117119
clc;
clearvars;
close all;
disp('Final Project (IE533)')
disp('Binary FSK - SIMULATION')

%% Propriedades de Modulacao
symbols_set = [0 1];

fc0 = 440;
fc1 = 4*fc0;
fs = 4*fc1;

%% Codigo corretor de erros
trellis = poly2trellis(3, [5 7]);

%% Sincronismo de quadros (MLS)
k = 10;

%% Fonte
% 1) Entrada de texto pelo usuario
% disp('Enter text:');
% input_str = input('','s');

% 2) Entrada de teste
input_str = 'abcdefghijklmnopqrstuvxzwyABCDEFGHIJKLMNOPQRSTUVXZWY1234567890';

%%
% Inicia um stopwatch timer
tic;
%% Vetor binario
%input_bin = str_source(input_str(i));
input_bin = str_source(input_str);

%% Mapeador
a = mapper(input_bin, symbols_set);
num_bits = length(a);
disp(['Number of bits sent:' num2str(num_bits)]);
%% Modulador
[s, t] = mod_fsk(a, fc0, fc1, fs, trellis, k);

figure()
plot(t, s)
title('s(t)')

%% Canal
% r = s;
% r = s + n;

% atraso no tempo
% ch_delay1 = round(100*rand(1)); % espacamento aleatorio
% ch_delay2 = round(100*rand(1));
% r = [zeros(1, ch_delay1) s zeros(1, ch_delay2)];

% pause(1); % delay

% Audio Record <---
numBits = 24;
numChannels = 1;
audio_rx = audiorecorder(fs, numBits, numChannels);
record(audio_rx);
pause(0.1); % delay

% B2B
sound(s, fs)

% input('Press ENTER to continue')
pause(t(end) + 0.1); % delay
stop(audio_rx)

r = getaudiodata(audio_rx);

figure()
plot(r)
title('r(t)')

%% Demodulador
z = demod_fsk(r, fc0, fc1, fs, trellis, k, length(a));
%% Demapeador
output_bin = demapper(z, symbols_set);

%% Destino
output_str = str_dest(output_bin);
% Mostra a saida de texto
disp('-----')
disp('Output text:');
fprintf('%s', output_str)

%%
disp(' ')
disp('-----')
% Calcula a taxa de erro de bit
% err = biterr(input_bin, output_bin);
err = biterr(a(1:length(z)), z);
disp(['Bit errors: ' num2str(err) ' (' num2str(100*err/num_bits) '%)'])
% Para o stopwatch timer
disp('Total time to compute:');
toc;