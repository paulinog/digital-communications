clc; clearvars; 
disp('FSK-2 example')
%% User parameters
numSymbol = 200;

fc0 = 440;
fc1 = 4*fc0;

fs = 4*fc1;

%% FEC
trellis = poly2trellis(3, [5 7]);
parity_ratio = 2;

%% Frame Sync
k = 8;
sync_bits = 2^k-1;
max_peak_pos = sync_bits * 0.8;

% TEST channel delays
ch_delay1 = 100*round(numSymbol*rand(1)); % random spacing
ch_delay2 = 100*round(numSymbol*rand(1));

%% Time vector
timestep = 1/fs;
T = 1/fc0;

frame_size = sync_bits + parity_ratio * numSymbol;

% MLS + FEC + Channel delays
% tmax = (ch_delay1 + sync_bits + parity_ratio * numSymbol + ch_delay2) * T;
% MLS + FEC
tmax = frame_size * T;
% FEC
% tmax = (parity_ratio * numSymbol) * T;

t = 0:timestep:tmax-timestep;

%% TX
a = randsrc(1, numSymbol, [0 1]);
a_enc = convenc(a, trellis);

sync_vec = double(mls(k, 1) > 0);
a_sync = [sync_vec a_enc]; % concatenate

h0 = sin(2*pi*fc0*t); % bit 0
h1 = sin(2*pi*fc1*t); % bit 1

w = heaviside(t) - heaviside(t - T);

m = upsample(a_sync, T*fs);
am = conv(m, w);
% plot(am(1:length(t)))

s = am(1:length(h0)) .* h1 ...
    + (1 - am(1:length(h1))) .* h0;
plot(t, s)

% sound(s, fs)

%% Channel 
% r = [zeros(1, ch_delay1) s zeros(1, ch_delay2)];
% r = [zeros(1, ch_delay1) s];
r = s;

frame_size_rx = frame_size;
tmax_rx = frame_size_rx * T;
t_rx = 0:timestep:tmax-timestep;

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

%% RX zero crossing counter
len_r = length(r);
c = zeros(1, len_r);
for i = 1 : len_r
    count = numel(find(t(zx) > (i-1)*T & t(zx) <= i*T));
    c(i) = count;
end
plot(t(1:numSymbol), c(1:numSymbol))
ylabel('c')

y_sync = double(c >= 5);
% stem(y_sync)

%% RX removing MLS
sync_vec2 = double(mls(k, 1) > 0);

self_corr = xcorr(y_sync, sync_vec2);
figure()
plot(self_corr(length(y_sync): length(y_sync) + frame_size));
% 
% start_frame = find(self_corr(length(y_sync): end) > max_peak_pos);
% 
% if (start_frame(1) + parity_ratio*numSymbol) > length(y_sync)
%     end_frame = length(y_sync);
% else
%     end_frame = (start_frame(1) + parity_ratio*numSymbol);
% end
% 
% y_enc = y_sync(start_frame(1) : end_frame);

%%

start_frame = 256;
% start_frame = ch_delay1 + (sync_bits+1);
end_frame = (start_frame-1) + (frame_size - sync_bits);

y_enc = y_sync(start_frame : end_frame);

%% FEC
y = vitdec(y_enc, trellis, 20,  'trunc', 'hard');

%% BER
err = sum(a ~= y);
if (err >0 )
    error([num2str(err) ' errors'])
else
    disp('no errors')
end

