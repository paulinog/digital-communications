%% Run Simulation
% Disciplina IE533A, 2s2019 (FEEC/UNICAMP)
% Guilherme Paulino, RA 117119
clc;
clearvars;
% close all;
disp('Final Project (IE533)')
%%  Parametros de Usuario
VERBOSE = false; % imprimir textos de log no Command Window a cada procedimento
OPEN_PLOT = false; % abrir janelas de graficos

%% Propriedades de Modulacao
% A) PAM-4
symbols_set = [-3, -1, 1, 3];

%% Fonte
% Entrada de texto pelo usuario
% disp('Enter text:');
% input_str = input('','s');

% Entrada de teste
% input_str = 'abcdefghijklmnopqrstuvxzwyABCDEFGHIJKLMNOPQRSTUVXZWY1234567890';
input_str = 'K';

% Inicia um stopwatch timer
tic;

%%
% example 16-QAM
% https://www.mathworks.com/help/comm/ug/pulse-shaping-using-a-raised-cosine-filter.html

err = 0;
for i = 1 : length(input_str)
    % Vetor binario
    input_bin = str_source(input_str(i));
    %% Mapeador
    a = mapper(input_bin, symbols_set);
    %% Modulador
    [s, t] = mod_pam4(a);
    %% Canal
    % r = s + n;
    %% Demodulador
    % z = demod_pam4(r);
    %% Demapeador
    output_bin = demapper(a, symbols_set);
    %% Destino
    err = err + biterr(input_bin, output_bin);
    output_str = str_dest(output_bin);
    % Mostra a saida de texto
    %disp('-----')
    %disp('Output text:');
    %disp(text_output)
    fprintf('%s', output_str)
end % end for-loop

%%
disp(' ')
disp('-----')
% Calcula a taxa de erro de bit
disp(['Bit error rate: ' num2str(err)])
% Para o stopwatch timer
disp('Total time to compute:');
toc;