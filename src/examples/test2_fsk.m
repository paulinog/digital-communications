clc; clearvars; 
disp('FSK-2 example')
%% User parameters
numSymbol = 10;

fc0 = 880;
fc1 = 2*fc0;

fs = 20* fc1;

%%
timestep = 1/fs;
T = 1/fc0;
tmax = numSymbol * T;
t = 0:timestep:tmax-timestep;
%% TX
a = randsrc(1, numSymbol, [0 1])

h0 = sin(2*pi*fc0*t);
h1 = sin(2*pi*fc1*t);

w = heaviside(t) - heaviside(t - T);

m = upsample(a, T*fs);
am = conv(m, w);

% plot(am)
%% TX sent
s = am(1:length(t)) .* h0 ...
    + (1 - am(1:length(t))) .* h1;

% figure();
% plot(t, s)

% sound(s, fs)

%% Channel 
% r = [zeros(1,100) mls s];
r = s;

%% RX
zci = @(v) find(v(:).*circshift(v(:), [-1 0]) <= 0);

zx = zci(r);

if (false)
    figure()
    plot(t, r, '-r')
    hold on
    plot(t(zx), r(zx), 'bp')
    hold off
    grid
    legend('Signal', 'Approximate Zero-Crossings')
end

c = zeros(1, numSymbol);
c_aux = 1;
for i = 1:numSymbol
    count = numel(find(t(zx) <= i*T));
    c(i) = count;
    c(i) = c(i) - c_aux;
    c_aux = count;
end

y = double(c <= 3)

