clc; clearvars; 
disp('FSK-2 example')
%% User parameters
numSymbol = 1000;

fc0 = 440;
fc1 = 2*fc0;

fs = 10*fc1;

%% FEC
trellis = poly2trellis(3, [5 7]);
parity_ratio = 2;

%% Frame Sync
k = 10;
sync_bits = 2^k-1;
max_peak_pos = sync_bits * 0.9;

% TEST channel delays
% ch_delay1 = round(numSymbol*rand(1));
% ch_delay2 = round(numSymbol*rand(1));

%% Time vector
timestep = 1/fs;
T = 1/fc0;
% tmax = (ch_delay1 + sync_bits + parity_ratio * numSymbol + ch_delay2) * T;
tmax = (sync_bits + parity_ratio * numSymbol) * T;
t = 0:timestep:tmax-timestep;

%% TX
a = randsrc(1, numSymbol, [0 1]);
a_enc = convenc(a, trellis);
a_sync = [double(mls(k, 1) > 0) a_enc];

h0 = sin(2*pi*fc0*t);
h1 = sin(2*pi*fc1*t);

w = heaviside(t) - heaviside(t - T);

m = upsample(a_sync, T*fs);
am = conv(m, w);
% plot(am)

s = am(1:length(t)) .* h0 ...
    + (1 - am(1:length(t))) .* h1;
% plot(t, s)

% sound(s, fs)

%% Channel 
% r = [zeros(1, ch_delay1) s zeros(1, ch_delay2)];
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

len_r = length(r);

% c = zeros(1, numSymbol * parity_ratio);
c = zeros(1, len_r);
c_aux = 1;
% for i = 1 : (numSymbol * parity_ratio)
for i = 1 : len_r
    count = numel(find(t(zx) <= i*T));
    c(i) = count;
    c(i) = c(i) - c_aux;
    c_aux = count;
end
c;
y_enc = double(c <= 3);

self_corr = xcorr(y_enc, double(mls(k, 1) > 0));
figure()
plot(self_corr)

start_frame = find(self_corr > max_peak_pos);
% vec_rx_i = vec_rx(start_frame(i) + 1 : start_frame(i) + tam);

% y = vitdec(y_enc, trellis, 20,  'trunc', 'hard');
% 
% %% BER
% err = sum(a ~= y);
% if (err >0 )
%     error([num2str(err) ' errors'])
% else
%     disp('no errors')
% end

