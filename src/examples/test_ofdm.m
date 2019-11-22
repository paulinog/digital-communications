close all;
clc; clearvars; 
disp('OFDM example')
%% User parameters
numSymbol = 16; % TODO: verificar se numSymbol for multiplo de 8

N = 8; % N-point FFT

% fc = 128; % frequencia da portadora
fc = 256;

fs = 8192; % frequencia de amostragem
% fs = 1024;

T = 1; % periodo do simbolo, em segundos
% R = N/T;

%% Time vector
timestep = T/(N*fs);
tmax = numSymbol*T/N;
t = 0:timestep:tmax-timestep;
t_pt = -tmax/2:timestep:tmax/2-timestep;

%% TX
a = randsrc(1, numSymbol, [-1 1]);

figure()
stem(a)
title('TX')
xlabel('samples')

% serial para paralelo
an = reshape(a, [N length(a)/N]);

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

%% BB
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

%%
figure()
hold on
plot(t, real(st))
plot(t, imag(st))
legend('real','imag')
xlabel('time')
ylabel('s(t)')


%% RF
% SRF = st.*exp(1j*2*pi*fc*t);

SRF = real(st).*cos(2*pi*fc*t) - imag(st).*sin(2*pi*fc*t);

figure()
hold on
plot(t, real(SRF))
plot(t, imag(SRF))
legend('real','imag')
xlabel('time')
ylabel('SRF')

%%
RRF_I = SRF .*cos(2*pi*fc*t);
RRF_Q = SRF .* -sin(2*pi*fc*t);

figure()
hold on
plot(t, RRF_I)
plot(t, RRF_Q)
legend('real','imag')
xlabel('time')
ylabel('RRF')

%%
LP_I = lowpass(RRF_I, fc/2, fs);
LP_Q = lowpass(RRF_Q, fc/2, fs);

figure()
hold on
plot(t, LP_I)
plot(t, LP_Q)
legend('real','imag')
xlabel('time')
ylabel('LP')

%%
% s(t) = a + jb = |s| * e^j<s = |s| * (cos(<s) + jsin(<s))
% s(t)*e^jx = (Re(s)*cos(x) - Im(s)*sin(x))
%             + j (Re(s)*sin(x) + Im(s)*cos(x))
%
% txbp = real(SRF)*cos(2*pi*fc*t) - imag(SRF)*sin(2*pi*fc*t);

%% Channel 
% noiseless
% r = st;

% real valued noise
% N0 = 0.8*max(st);
% n = N0*rand(1, length(t));
% r = st + n;

% complex noise
% snr = -40; %dB
snr = 0;
powerDB = 10*log10(var(st));
noiseVar = 10.^(0.1*(powerDB-snr)); 
r = awgn(st, snr);

figure()
subplot(211)
plot(t,abs(r))
xlabel('time')
ylabel('|r + n|')
subplot(212)
plot(t,angle(r))
xlabel('time')
ylabel('\angle{r + n}')


%% RX OFDM
% rxbb=hilbert(SRF).*exp(-j2pift)

tic
rt = conv(r, p(t_pt), 'same');
disp('Time to compute convolution:')
toc

figure()
subplot(211)
plot(t,abs(rt))
xlabel('time')
ylabel('|r(t)|')
subplot(212)
plot(t,angle(rt))
xlabel('time')
ylabel('\angle{r(t)}')

% rk_up = rt(t = kT);
rk = downsample(rt, fs);

% serial para paralelo
rkn = reshape(rk, [N length(rk)/N]);

yn = fft( rkn, N );
% y_p = sign(real(fft( rkn, N )));

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
    error(['Total of errors: ' num2str(err) ' (' ...
           num2str(100*err/numSymbol) '%)']);
else
    disp('no errors');
end

