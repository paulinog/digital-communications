function [y] = demod_fsk(r, fc0, fc1, fs, trellis, k, numSymbol)

timestep = 1/fs;
T = 1/fc0;

%% Frame Sync
sync_bits = 2^k-1;

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
if (true)
    figure()
    plot(t_rx, r, '-r')
    hold on
    plot(t_rx(zx), r(zx), 'bp')
    hold off
    grid
    title('Zero-crossings counter')
    legend('Signal', 'Approximate Zero-Crossings')
    ylabel('r(t)')
    xlabel('time')
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
figure()
plot(t_rx, c)
ylabel('c')

% disp(['max c= ' num2str(max(c))])
% disp(['min c= ' num2str(min(c))])

y_sync = double(c >= 4);

% figure()
% stem(y_sync)

%% RX removing MLS
sync_vec = double(mls(k, 1) > 0);

self_corr = xcorr(y_sync, sync_vec);
figure()
plot(self_corr);
title('Cross correlation')
ylabel('Rs')
xlabel('time')

max_peak_pos = max(self_corr);
if(length(max_peak_pos) ~= 1)
    disp('max_peak_pos = ')
    disp(max_peak_pos)
end

start_frame = find(self_corr(length(y_sync): end) == max_peak_pos) + sync_bits;

if(isempty(start_frame) || (length(start_frame) > 1))
    disp('start_frame =')
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

%% FEC
y = vitdec(y_enc, trellis, 20,  'trunc', 'hard');

end