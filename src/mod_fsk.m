function [s,t] = mod_fsk(a, fc0, fc1, fs, trellis, k)
numSymbol = length(a);

%% FEC

parity_ratio = 2;

%% Frame Sync
sync_bits = 2^k-1;

%% Time vector
timestep = 1/fs;
T = 1/fc0;

frame_size = sync_bits + parity_ratio * numSymbol;

% MLS + FEC
tmax = frame_size * T;

t = 0:timestep:tmax-timestep;

%% TX
a_enc = convenc(a, trellis);

sync_vec = double(mls(k, 1) > 0);
a_sync = [sync_vec a_enc]; % concatenate

h0 = sin(2*pi*fc0*t); % bit 0
h1 = sin(2*pi*fc1*t); % bit 1

w = heaviside(t) - heaviside(t - T); % window function

m = upsample(a_sync, T*fs);
am = conv(m, w);
% figure()
% plot(am(1:length(t)))

s = am(1:length(h1)) .* h1 ...
    + (1 - am(1:length(h0))) .* h0;
% figure()
% plot(t, s)
end
