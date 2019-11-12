%% Run Transmitter
% Disciplina IE533A, 2s2019 (FEEC/UNICAMP)
% Guilherme Paulino, RA 117119
clc;
clearvars;
% close all;
disp('Final Project (IE533)')
disp('Binary FSK - TRANSMITTER')

%% Propriedades de Modulacao
symbols_set = [0 1];

fc0 = 440;
fc1 = 4*fc0;
fs = 4*fc1;

%% Codigo corretor de erros
trellis = poly2trellis(3, [5 7]);

%% Sincronismo de quadros (MLS)
k = 8;

%% Fonte
% Entrada de texto pelo usuario
disp('Enter text:');
input_str = input('','s');

%%
% Inicia um stopwatch timer
tic;
%% Vetor binario
%input_bin = str_source(input_str(i));
input_bin = str_source(input_str);

%% Mapeador
a = mapper(input_bin, symbols_set);
disp(['Number of bits sent:' num2str(length(a))]);
%% Modulador
s = mod_fsk(a, fc0, trellis, k);

%% Transmitir
sound(s, fs)

%% Para o stopwatch timer
disp('-----')
disp('Total time to send:');
toc;