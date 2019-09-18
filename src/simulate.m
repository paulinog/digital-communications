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
% Text input by the user
disp('Enter text:');
array_bin = str_source();
% Start a stopwatch timer
tic;

%% Mapper
a = mapper(array_bin, symbols_set);

%% Modulator
% s = mod_pam4(a);
%% Channel
% r = s+n;
z = a;
%% Demodulator
%% Demapper
out_array = demapper(z, symbols_set);

%% Destination
text_output = str_dest(out_array);
% Stop stopwatch timer
disp('-----')
disp('Total time to compute:');
toc;
% Displays text output
disp('-----')
disp('Output text:');
disp(text_output)
