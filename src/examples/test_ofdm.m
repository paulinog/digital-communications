close all;
clc; clearvars; 
disp('OFDM example')
%% User parameters
numSymbol = 32; % TODO: verificar se numSymbol for multiplo de 8

N = 8; % N-point FFT

% fc = 128; % frequencia da portadora
% fc = 256;
fc = 1e3;

fs = 8192; % frequencia de amostragem
% fs = 1024;

T = 1; % periodo do simbolo, em segundos
% R = N/T;

%% Parametros do FEC
enable_FEC = true;
trellis = poly2trellis(3, [5 7]);
parity_ratio = 2;

%% Parametros do MLS
enable_MLS = true;
k = 8;
sync_bits = 2^k-1;

% TEST channel delays
ch_delay1 = 100*round(100*numSymbol*rand(1)); % random spacing
ch_delay2 = 100*round(100*numSymbol*rand(1));

%% Vetor Tempo
timestep = T/(N*fs);

frame_size = numSymbol;
if enable_FEC
    frame_size = frame_size * parity_ratio;
end
if enable_MLS
    frame_size = frame_size + 2*sync_bits;
end
num_padding = zero_padding(frame_size, N);
if num_padding > 0
    frame_size = frame_size + num_padding;
end
tmax = frame_size*T/N;
t = 0:timestep:tmax-timestep;
t_pt = -tmax/2:timestep:tmax/2-timestep;

%% TX
a = randsrc(1, numSymbol, [0 1]);

figure()
stem(a)
title('TX')
xlabel('samples')

%% TX FEC
if enable_FEC
    a_enc = convenc(a, trellis);
else
    a_enc = a;
end

%% TX MLS
if 0% enable_MLS
    sync_vec = double(mls(k, 1) > 0);
%     a_sync = [sync_vec a_enc]; % concatenate
    a_sync = [sync_vec a_enc sync_vec]; % concatenate
else
    a_sync = a_enc;
end

%% TX BPSK
a_mod = (a_sync >= 0.5) - (a_sync < 0.5);

%% TX Serial para paralelo
num_padding = zero_padding(length(a_mod), N);

if num_padding > 0
    a_mod = [a_mod zeros(1, num_padding)];
end

an = reshape(a_mod, [N length(a_mod)/N]);

%% TX OFDM
skn = ifft(an, N);

figure()
stem(abs(skn)')
title('IFFT')
ylabel('|s_k|')
xlabel(['samples / ' num2str(N)])
figure()
stem(angle(skn)')
title('IFFT')
ylabel('\angle{s_k}')

%% TX Paralelo para serial
sk = reshape(skn, [1 size(skn, 1)*size(skn, 2)]);

% Add by Marcos
if  enable_MLS
    sync_vec = double(mls(k, 1) > 0);
%     a_sync = [sync_vec a_enc]; % concatenate
    sk_sync = [sync_vec sk sync_vec]; % concatenate
else
    sk_sync = sk;
end

%% TX Upsample
% sk_up = upsample(sk, fs);
sk_up = upsample(sk_sync, fs); %Add by Marcos


figure()
plot(t, abs(sk_up))
title('S_k upsampled')
xlabel('time (in seconds)')

% p = @(t) sqrt(T/N) * sin(pi*N*t/T) ./ (pi*t);
p = @(t) sqrt(N/T) * sinc(N*t/T);

figure();
pt = p(t_pt);
plot(t_pt,pt);
xlabel('time')
ylabel('p(t)')

%% TX Base Band
st = conv (sk_up, pt, 'same');

figure()
subplot(211)
plot(t,abs(st))
xlabel('time')
ylabel('|s(t)|')
subplot(212)
plot(t,angle(st))
xlabel('time')
ylabel('\angle{s(t)}')

%% TX sinal complexo
figure()
hold on
plot(t, real(st))
plot(t, imag(st))
legend('real','imag')
xlabel('time')
ylabel('s(t)')

%% TX Pass Band
SRF = real(st).*cos(2*pi*fc*t) - imag(st).*sin(2*pi*fc*t);

figure()
hold on
plot(t, real(SRF))
plot(t, imag(SRF))
legend('real','imag')
xlabel('time')
ylabel('SRF')

%% Channel 
close all

% noiseless
r = SRF;
% r = [zeros(1, ch_delay1) SRF];
% r = [zeros(1, ch_delay1) SRF zeros(1, ch_delay2)];
len_r = length(r);
%%

t_rx = 0:timestep:(len_r-1)*timestep;

%%
% real valued noise
% N0 = 0.8*max(st);
% n = N0*rand(1, length(t));
% r = SRF + n;

% complex noise
% snr = -40; %dB
% snr = 0;
% powerDB = 10*log10(var(SRF));
% noiseVar = 10.^(0.1*(powerDB-snr)); 
% r = awgn(SRF, snr);
%%
figure()
subplot(211)
plot(t_rx,abs(r))
xlabel('time')
ylabel('|r + n|')
subplot(212)
plot(t_rx,angle(r))
xlabel('time')
ylabel('\angle{r + n}')

%% RX RF
phaseRF = 0;
RRF_I = r .*cos(2*pi*fc*t_rx + phaseRF);
RRF_Q = r .* -sin(2*pi*fc*t_rx + phaseRF);

figure()
hold on
plot(t_rx, RRF_I)
plot(t_rx, RRF_Q)
legend('real','imag')
xlabel('time')
ylabel('RRF')

%% RX LP
[num, den] = butter(10, fc*2*timestep, 'low');
LP_I = filtfilt(num, den, RRF_I) * 2;
LP_Q = filtfilt(num, den, RRF_Q) * 2;

figure()
hold on
plot(t_rx, LP_I)
plot(t_rx, LP_Q)
legend('real','imag')
xlabel('time')
ylabel('LP')

%% RX Complex signal
rt = LP_I + 1j*LP_Q;

figure()
subplot(211)
plot(t_rx,abs(rt))
xlabel('time')
ylabel('|r(t)|')
subplot(212)
plot(t_rx,angle(rt))
xlabel('time')
ylabel('\angle{r(t)}')

%% RX OFDM
% rk_up = rt(t = kT);
rk = downsample(rt, fs);

num_padding = zero_padding(length(rk), N);

if num_padding > 0
    rk = [rk zeros(1, num_padding)];
end

% serial para paralelo
rkn = reshape(rk, [N length(rk)/N]);

% y_p = sign(real(fft( rkn, N )));
yn = fft( rkn, N );

%% RX Paralelo para serial
y = reshape(yn, [1 size(yn, 1)*size(yn, 2)]);

%% RX Slicer
z_sync = (y > 0);

figure()
stem(z_sync)
title('RX')
xlabel('samples')
ylabel('z')

%% RX removing MLS
if enable_MLS
    sync_vec2 = double(mls(k, 1) > 0);
    self_corr = xcorr(z_sync, sync_vec2);
    figure()
    plot(self_corr);
    xlim([length(z_sync)-length(sync_vec) length(z_sync)+frame_size])
    title('Cross correlation')
    ylabel('R')
    xlabel('sample')

    max_peak_pos = max(self_corr);
%     if(length(max_peak_pos) ~= 1)
%         warning('max_peak_pos = ')
%         disp(max_peak_pos)
%     end
%     start_frame = find(self_corr(length(z_sync): end) == max_peak_pos) + sync_bits;
    start_frame = min(find(self_corr(length(z_sync): end) == max_peak_pos)) + sync_bits;
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
%     end_frame = (start_frame-1) + (frame_size - sync_bits - num_padding);
    end_frame = max(find(self_corr(length(z_sync): end) == max_peak_pos)) - 1;
    z_enc = z_sync(start_frame : end_frame);
else
    z_enc = z_sync(1:frame_size);
end

%% RX FEC
if enable_FEC
    z = vitdec(z_enc, trellis, 20,  'trunc', 'hard');
else
    z = z_enc;
end

%% RX BER
err = sum(a ~= z);
if (err >0 )
    error(['Total of errors: ' num2str(err) ' (' ...
           num2str(100*err/frame_size) '%)']);
else
    disp('no errors');
end

