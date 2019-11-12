function [y] = demod_fsk(r, numSymbol, fc0, trellis, k)
fc1 = 4*fc0;
fs = 4*fc1;

timestep = 1/fs;
T = 1/fc0;

%% Frame Sync
sync_bits = 2^k-1;
max_peak_pos = 120;

%% FEC
parity_ratio = 2;

%% Time Vector
frame_size = sync_bits + parity_ratio * numSymbol;
len_r = length(r);
tmax_rx = len_r/fs;
t_rx = 0 : timestep : tmax_rx - timestep;
% plot(t_rx, r)

%% RX
zci = @(v) find(v(:).*circshift(v(:), [-1 0]) <= 0);
zx = zci(r);
if (false)
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
    if c(i) > 10
        c(i) = 0;
    end
end
% figure()
% plot(t_rx, c)
% ylabel('c')

y_sync = double(c >= 5);

% figure()
% stem(y_sync)

%% RX removing MLS
sync_vec2 = double(mls(k, 1) > 0);

self_corr = xcorr(y_sync, sync_vec2);
% figure()
% plot(self_corr);

start_frame = find(self_corr(length(y_sync): end) > max_peak_pos) + sync_bits;
% 
% if (start_frame(1) + parity_ratio*numSymbol) > length(y_sync)
%     end_frame = length(y_sync);
% else
%     end_frame = (start_frame(1) + parity_ratio*numSymbol);
% end
% 
% y_enc = y_sync(start_frame(1) : end_frame);

end_frame = (start_frame-1) + (frame_size - sync_bits);
y_enc = y_sync(start_frame : end_frame);

%% FEC
y = vitdec(y_enc, trellis, 20,  'trunc', 'hard');

end