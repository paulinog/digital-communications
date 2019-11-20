%% MLS Test
clc; clearvars;
disp ('MLS Test')

k = 8;
tam = 1000;

num_bits = 2^k-1; % number of bits

%% TX
disp(['number of bits sequence = ' num2str(num_bits)])
sync_vec_tx = mls(k, 1);

vec1 = randsrc(1, tam, [1 -1]);
vec2 = randsrc(1, tam, [1 -1]);
vec3 = randsrc(1, tam, [1 -1]);

vec_tx = [vec1 sync_vec_tx vec2 sync_vec_tx vec3];

%% Channel
vec_rx = vec_tx;

%% RX

sync_vec_rx = mls(k, 1);

self_corr = xcorr(sync_vec_rx, vec_rx);
plot(self_corr)

max_peak_pos = max(self_corr)-1;
start_frame = find(self_corr > max_peak_pos);

for i = 1 : length(start_frame)
    vec_rx_i = vec_rx(start_frame(i) + 1 : start_frame(i) + tam);

    % TEST RX
    if (vec2 == vec_rx_i)
        disp([num2str(i) ': equal to vec2'])
    end
    if (vec3 == vec_rx_i)
        disp([num2str(i) ': equal to vec3'])
    end
end
