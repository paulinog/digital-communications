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

fc0 = 110;
fc1 = 4*fc0;
fs = 8*fc1;

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

%% Inicia um stopwatch timer
tic;

%% Vetor binario
input_bin = str_source(input_str);

%% Mapeador
a = mapper(input_bin, symbols_set);
num_bits = length(a);
disp(['Number of bits sent:' num2str(num_bits)]);
%% Modulador
[s, t] = mod_fsk(a, fc0, fc1, fs, trellis, k);

figure()
plot(t, s)
title('TX')
ylabel('s(t)')
xlabel('time')

%% Canal
% r = s + n;
% r=s;

% grava o audio
numBits = 16;
numChannels = 1;
audio_rx = audiorecorder(fs, numBits, numChannels);
record(audio_rx);
pause(0.1); % delay

% transmit ->
sound(s, fs)

pause(t(end) + 0.1); % delay
stop(audio_rx)

% receive <-
r = getaudiodata(audio_rx)';

figure()
t_ch = 0 : 1/fs : (length(r)-1)/fs ;
plot(t_ch, r)
title('RX')
ylabel('r(t)')
xlabel('time')

%% RX Demodulador
z = demod_fsk(r, fc0, fc1, fs, trellis, k, length(a));
%% Demapeador
output_bin = demapper(z, symbols_set);

%% RX Destino
output_str = str_dest(output_bin);
% Mostra a saida de texto
disp('-----')
disp('Output text:');
fprintf('%s', output_str)

%% BER
disp(' ')
disp('-----')
% Calcula a taxa de erro de bit
if (length(a) == length(z))
    err = biterr(a, z);
    disp(['Bit errors: ' num2str(err) ' (' num2str(100*err/num_bits) '%)'])
else
    disp('Bit errors: z has a different size')
end
% Para o stopwatch timer
disp('Total time to compute:');
toc;