clc; clearvars; 
disp('FSK-2 example')
%% User parameters
numSymbol = 1000;

fc0 = 440;
fc1 = 2*fc0;

fs = 100* fc0;

%%
timestep = 1/fs;
T = 1/fc0;
tmax = numSymbol * T;
t = 0:timestep:tmax-timestep;
%% TX
a = randsrc(1, numSymbol, [0 1]);

h0 = sin(2*pi*fc0*t);
h1 = sin(2*pi*fc1*t);

w = heaviside(t) - heaviside(t - T);

m = upsample(a, T*fs);
am = conv(m, w);

% plot(am)
%% RX
s = am(1:length(t)) .* h0 ...
    + (1 - am(1:length(t))) .* h1;
plot(t, s)

sound(s, fs)