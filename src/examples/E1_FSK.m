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
max_peak_pos = (sync_bits) * 0.9;

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

w = heaviside(t) - heaviside(t - T); % window function

m = upsample(a_sync, T*fs);
am = conv(m, w);
% plot(am(1:length(t)))

s = am(1:length(h1)) .* h1 ...
    + (1 - am(1:length(h0))) .* h0;
% plot(t, s)

% input('Press ENTER to continue')
% sound(s, fs)


%% Audio Record
% FS = 8192;
numBits = 8;
numChannels = 1;
% audioRecordTime = 10; %in seconds
% audioWaitUnit = 0.1; %in seconds

audio_rx = audiorecorder(8*fs, numBits, numChannels);

input('Press ENTER to continue')

disp('starting of record')
record(audio_rx);

pause(0.5);
sound(s, fs)

%pause(audioRecordTime + 1);

input('Press ENTER to continue')

disp('endding record')
stop(audio_rx)


%% Channel 
% r = [zeros(1, ch_delay1) s zeros(1, ch_delay2)];
% r = [zeros(1, ch_delay1) s ];
% r = s;
r = getaudiodata(audio_rx);

plot(r)

len_r = length(r);
% frame_size_rx = len_r/(fs*T);
tmax_rx = len_r/fs;
t_rx = 0 : timestep : tmax_rx - timestep;
% plot(t_rx, r)

%% RX
% zci = @(v) find(v(:).*circshift(v(:), [-1 0]) <= 0);
zci = @(v) find(v(:).*circshift(v(:), 0) <= 0);

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

%% RX zero crossing counter

c = zeros(1, len_r);
for i = 1 : len_r
    count = numel(find(t_rx(zx) > (i-1)*T & t_rx(zx) <= i*T));
    c(i) = count;
%     if c(i) >= 8
%         c(i) = 0;
%     end
end
plot(t_rx, c)
% plot(t(1:frame_size), c(1:frame_size))
% ylabel('c')

% y_sync = double(c >= 5);

% figure()
% stem(y_sync)

%% RX removing MLS
sync_vec2 = double(mls(k, 1) > 0);

self_corr = xcorr(y_sync, sync_vec2);
% figure()
% plot(self_corr(length(y_sync): length(y_sync) + frame_size));
% plot(self_corr);
% 

% start_frame = find(self_corr(length(y_sync): end) > max_peak_pos) + sync_bits
start_frame = find(self_corr(length(y_sync): end) > 100) + sync_bits

% 
% if (start_frame(1) + parity_ratio*numSymbol) > length(y_sync)
%     end_frame = length(y_sync);
% else
%     end_frame = (start_frame(1) + parity_ratio*numSymbol);
% end
% 
% y_enc = y_sync(start_frame(1) : end_frame);

%%

% start_frame = (sync_bits+1);
% start_frame = (sync_bits+1) + ch_delay1;

%start_frame = sync_bits + 875;
% start_frame = 1020 + 255;


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

