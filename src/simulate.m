%% Run Simulation
% Disciplina IE533A, 2s2019 (FEEC/UNICAMP)
% Guilherme Paulino, RA 117119
% clc;
clearvars;
close all;
disp('Final Project (IE533)')
%%  Parametros de Usuario
VERBOSE = false; % imprimir textos de log no Command Window a cada procedimento
OPEN_PLOT = false; % abrir janelas de gráficos

%% Propriedades de Modulacao
% A) PAM-4
symbols_set = [-3, -1, 1, 3];

% while true
%% Fonte
% Entrada de texto pelo usuario
disp('Enter text:');
array_bin = str_source();
% Inicia um stopwatch timer
tic;

%% Mapeador
a = mapper(array_bin, symbols_set);
%% Filtro Formatador de Pulso
% https://www.mathworks.com/help/comm/ug/pulse-shaping-using-a-raised-cosine-filter.html
%% Moduador
% s = mod_pam4(a);
%% Canal
% r = s+n;
z = a;
%% Demodulador
% z = demod_pam4(r);
%% Demapeador
out_array = demapper(z, symbols_set);

%% Destino
text_output = str_dest(out_array);
% Para o stopwatch timer
disp('-----')
disp('Total time to compute:');
toc;
% Mostra a saida de texto
disp('-----')
disp('Output text:');
disp(text_output)
% end % end-while