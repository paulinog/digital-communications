close all;
clc; clearvars; 
disp('OFDM example')
%% User parameters
numSymbol = 80; % TODO: verificar se numSymbol for multiplo de 8

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

%% Vetor Tempo
timestep = T/(N*fs);

frame_size = numSymbol;
if enable_FEC
    frame_size = frame_size * parity_ratio;
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

%% TX BPSK
a_mod = (a_enc >= 0.5) - (a_enc < 0.5);

%% TX Serial para paralelo
an = reshape(a_mod, [N length(a_enc)/N]);

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

%% TX Upsample
sk_up = upsample(sk, fs);

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
% noiseless
r = SRF;

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

figure()
subplot(211)
plot(t,abs(r))
xlabel('time')
ylabel('|r + n|')
subplot(212)
plot(t,angle(r))
xlabel('time')
ylabel('\angle{r + n}')

%% RX RF
RRF_I = r .*cos(2*pi*fc*t);
RRF_Q = r .* -sin(2*pi*fc*t);

figure()
hold on
plot(t, RRF_I)
plot(t, RRF_Q)
legend('real','imag')
xlabel('time')
ylabel('RRF')

%% RX LP
[num, den] = butter(10, fc*2*timestep);
LP_I = filtfilt(num, den, RRF_I) * 2;
LP_Q = filtfilt(num, den, RRF_Q) * 2;

figure()
hold on
plot(t, LP_I)
plot(t, LP_Q)
legend('real','imag')
xlabel('time')
ylabel('LP')

%% RX Complex signal
rt = LP_I + 1j*LP_Q;

figure()
subplot(211)
plot(t,abs(rt))
xlabel('time')
ylabel('|r(t)|')
subplot(212)
plot(t,angle(rt))
xlabel('time')
ylabel('\angle{r(t)}')

%% RX OFDM
% rk_up = rt(t = kT);
rk = downsample(rt, fs);

% serial para paralelo
rkn = reshape(rk, [N length(rk)/N]);

% y_p = sign(real(fft( rkn, N )));
yn = fft( rkn, N );

%% RX Paralelo para serial
y = reshape(yn, [1 size(yn, 1)*size(yn, 2)]);

%% RX Slicer
z_enc = (y > 0);

figure()
stem(z_enc)
title('RX')
xlabel('samples')
ylabel('z')

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

