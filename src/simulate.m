%% Start Simulation
% Disciplina IE533A/2s2019
% Guilherme Paulino, RA 117119
clc; clearvars; close all;
disp('Final Project (IE533)')
%%  User parameters
VERBOSE = false;

%% Modulation Proprierties (4-PAM)
symbols_set = [-3, -1, 1, 3];

%% Source
disp('Enter text:');
array_bin = str_source();

%% Mapper
a = mapper(array_bin, symbols_set);

%% Modulator
%% Channel
z = a;
%% Demodulator
%% Demapper
out_array = demapper(z, symbols_set);

%% Destination
text_output = str_dest(out_array);
disp('-----')
disp('Output text:');
disp(text_output)
