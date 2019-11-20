close all;
clc; clearvars; 
disp('OFDM example')
%% User parameters
numSymbol = 20;
fc0 = 440;
fc1 = 4*fc0;
fs = 4*fc1;

%% FEC Parameters
enable_FEC = true;

trellis = poly2trellis(3, [5 7]);
parity_ratio = 2;

%% MLS Parameters
enable_MLS = false;

k = 8;
sync_bits = 2^k-1;

% TEST channel delays
% ch_delay1 = 100*round(numSymbol*rand(1)); % random spacing
% ch_delay2 = 100*round(numSymbol*rand(1));

%% Time vector
timestep = 1/fs;
T = 1/fc0;

frame_size = numSymbol;
if enable_FEC
    frame_size = frame_size * parity_ratio;
end
if enable_MLS
    frame_size = frame_size + sync_bits;
end

% MLS + FEC + Channel delays
% tmax = (ch_delay1 + sync_bits + parity_ratio * numSymbol + ch_delay2) * T;
% MLS + FEC
tmax = frame_size * T;
% FEC
% tmax = (parity_ratio * numSymbol) * T;

t = 0:timestep:tmax-timestep;

%% TX
a = randsrc(1, numSymbol, [0 1]);

%% TX FEC
if enable_FEC
    a_enc = convenc(a, trellis);
else
    a_enc = a;
end

%% TX MLS
if enable_MLS
    sync_vec = double(mls(k, 1) > 0);
    a_sync = [sync_vec a_enc]; % concatenate
else
    a_sync = a_enc;
end

%% TX FSK-2
h0 = sin(2*pi*fc0*t); % bit 0
h1 = sin(2*pi*fc1*t); % bit 1
w = heaviside(t) - heaviside(t - T); % window function
m = upsample(a_sync, T*fs);
am = conv(m, w);
% plot(am(1:length(t)))
s = am(1:length(h1)) .* h1 ...
    + (1 - am(1:length(h0))) .* h0;
% plot(t, s)
% sound(s, fs)

%% TX BPSK

%% Channel 
r = s;
% r = [zeros(1, ch_delay1) s ];
% r = [zeros(1, ch_delay1) s zeros(1, ch_delay2)];

len_r = length(r);
% frame_size_rx = len_r/(fs*T);
tmax_rx = len_r/fs;
t_rx = 0 : timestep : tmax_rx - timestep;
plot(t_rx, r)

%% RX FSK-2
zci = @(v) find(v(:).*circshift(v(:), [-1 0]) <= 0);
zx = zci(r);
if (true)
    figure()
    plot(t_rx, r, '-r')
    hold on
    plot(t_rx(zx), r(zx), 'bp')
    hold off
    grid
    legend('Signal', 'Approximate Zero-Crossings')
end
c = zeros(1, len_r);
for i = 1 : len_r
    count = numel(find(t_rx(zx) > (i-1)*T & t_rx(zx) <= i*T));
    c(i) = count;
    if c(i) > 10
        c(i) = 0;
    end
end
% plot(t(1:frame_size), c(1:frame_size))
plot(t_rx, c)
xlim([0 frame_size/fs])
title('Zero-crossings counter per Period')
ylabel('counts / period')

y_sync = double(c >= 5);

figure()
stem(y_sync)

%% RX removing MLS
if enable_MLS
    sync_vec2 = double(mls(k, 1) > 0);
    self_corr = xcorr(y_sync, sync_vec2);
    figure()
    % plot(self_corr(length(y_sync): length(y_sync) + frame_size));
    plot(self_corr);
    xlim([length(y_sync)-length(sync_vec) length(y_sync)+frame_size])
    title('Cross correlation')
    ylabel('R')
    xlabel('sample')

    max_peak_pos = max(self_corr);
    if(length(max_peak_pos) ~= 1)
        warning('max_peak_pos = ')
        disp(max_peak_pos)
    end
    start_frame = find(self_corr(length(y_sync): end) == max_peak_pos) + sync_bits;
    if(isempty(start_frame) || (length(start_frame) > 1))
        warning('start_frame =')
        disp(start_frame)

        if isempty(start_frame)
            start_frame = 1;
        end
        if (length(start_frame) > 1)
            start_frame = start_frame(1);
        end
    end
    end_frame = (start_frame-1) + (frame_size - sync_bits);
    y_enc = y_sync(start_frame : end_frame);
else
    y_enc = y_sync(1:frame_size);
end
%% Decode FEC
if enable_FEC
    y = vitdec(y_enc, trellis, 20,  'trunc', 'hard');
else
    y = y_enc;
end
%% RX BER
err = sum(a ~= y);
if (err >0 )
    error([num2str(err) ' errors'])
else
    disp('no errors')
end

