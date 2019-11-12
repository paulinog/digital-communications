%% Run Simulation
% Disciplina IE533A, 2s2019 (FEEC/UNICAMP)
% Guilherme Paulino, RA 117119
% clc;
clearvars;
% close all;
disp('Final Project (IE533)')
disp('Binary FSK')
%%  Parametros de Usuario
VERBOSE = false; % imprimir textos de log no Command Window a cada procedimento
OPEN_PLOT = false; % abrir janelas de graficos

%% Propriedades de Modulacao

symbols_set = [0 1];

fc0 = 440;
% fc1 = 4*fc0;
% fs = 4*fc1;

%% Codigo corretor de erros
trellis = poly2trellis(3, [5 7]);

%% Sincronismo de quadros (MLS)
k = 8;

%% Fonte
% Entrada de texto pelo usuario
%disp('Enter text:');
% input_str = input('','s');

% Entrada de teste
input_str = 'abcdefghijklmnopqrstuvxzwyABCDEFGHIJKLMNOPQRSTUVXZWY1234567890';
% input_str = 'K';

% Inicia um stopwatch timer
tic;

err = 0;
% for i = 1 : length(input_str)
    % Vetor binario
    %input_bin = str_source(input_str(i));
    input_bin = str_source(input_str);
    %% Mapeador
    a = mapper(input_bin, symbols_set);
    numSymbol = length(a);
    %% Modulador
    [s, t] = mod_fsk(a, fc0, trellis, k);
    
    %% Canal
    %r = s;
    % r = s + n;
    
    % atraso no tempo
    ch_delay1 = round(numSymbol*rand(1)); % espacamento aleatorio
    ch_delay2 = round(numSymbol*rand(1));
    r = [zeros(1, ch_delay1) s zeros(1, ch_delay2)];
    
    %% Demodulador
    z = demod_fsk(r, numSymbol, fc0, trellis, k);
    %% Demapeador
    output_bin = demapper(z, symbols_set);
    
    %% Destino
    err = err + biterr(input_bin, output_bin);
    output_str = str_dest(output_bin);
    % Mostra a saida de texto
    %disp('-----')
    %disp('Output text:');
    %disp(text_output)
    fprintf('%s', output_str)
% end % end for-loop

%%
disp(' ')
disp('-----')
% Calcula a taxa de erro de bit
disp(['Bit errors: ' num2str(err)])
% Para o stopwatch timer
disp('Total time to compute:');
toc;