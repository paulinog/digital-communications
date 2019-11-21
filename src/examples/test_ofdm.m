close all;
clc; clearvars; 
disp('OFDM example')
%% User parameters
numSymbol = 80; % TODO: verificar se numSymbol for multiplo de 8

N = 8; % N-point FFT

fc = 128; % frequencia da portadora
fs = 8192; % frequencia de amostragem

T = 2/fc;
R = N/T;

%% Time vector
timestep = T/fs;
tmax = numSymbol*T;
t = 0:timestep:tmax-timestep;
t_pt = -tmax/2:timestep:tmax/2-timestep;

%% TX
a = randsrc(1, numSymbol, [-1 1]);

figure()
stem(a)
title('TX')
xlabel('samples')

% serial para paralelo
an = reshape(a, [8 length(a)/8]);

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

% paralelo para serial
sk = reshape(skn, [1 size(skn, 1)*size(skn, 2)]);

txbb = upsample(sk, fs);

figure()
plot(t, abs(txbb))
title('S_k upsampled')
xlabel('time (in seconds)')

% p = @(t) sqrt(T/N) * sin(pi*N*t/T) ./ (pi*t);
p = @(t) sqrt(N/T) * sinc(t/T);

figure();
pt = p(t_pt);
plot(t_pt,pt);
xlabel('time')
ylabel('p(t)')

%%
% st = conv (real(txbb), p(t), 'same');
st = conv (txbb, pt, 'same');

figure()
subplot(211)
plot(t,abs(st))
xlabel('time')
ylabel('|s(t)|')
subplot(212)
plot(t,angle(st))
xlabel('time')
ylabel('\angle{s(t)}')

%SRF = st.*exp(j*2*pi*f*t)
SRF = st.*exp(j*2*pi*fc*t);

figure()
hold on
plot(real(SRF))
plot(imag(SRF))

figure()
subplot(211)
plot(t,abs(SRF))
xlabel('time')
ylabel('|SRF|')
subplot(212)
plot(t,angle(SRF))
xlabel('time')
ylabel('\angle{SRF}')


%% Channel 
% rt = st;

% figure()
% hold on
% stem(real(rt))
% stem(imag(rt))
% hold off
% title('Complex Channel')
% ylabel('r(n)')
% xlabel('samples')
% legend('real', 'imaginary')

% plot(t_rx, r)
% ylabel('r(t)')
% xlabel('time')

%% RX OFDM
% rxbb=hilbert(SRF).*exp(-j2pift);
% rxbb = rt;
% rx = conv(rxbb, pt, 'same');
% y_p = sign(real(fft( rk, numSC )));

rkn = skn; % bypass modulation
yn = fft( rkn, N );

% paralelo para serial
y = reshape(yn, [1 size(yn, 1)*size(yn, 2)]);

% slicer
z = (y >= 0) - (y < 0);

figure()
stem(z)
title('RX')
xlabel('samples')
ylabel('z')

%% RX BER
err = sum(a ~= z);
if (err >0 )
    error([num2str(err) ' errors'])
else
    disp('no errors')
end

