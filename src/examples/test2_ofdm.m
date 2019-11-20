close all;
clc; clearvars; 
disp('OFDM example')
%% User parameters
numSymbol = 800;
numSC = 8;

fc = 430;
fs = 8192;

% TODO: verificar se numSymbol for multiplo de 8

%% Time vector
timestep = 1/fs;
T = 1/fc;

frame_size = numSymbol;
tmax = frame_size * T;

t = 0:timestep:tmax-timestep;

%% TX
a = randsrc(1, numSymbol, [-1 1]);

figure()
stem(a)
title('TX')
xlabel('samples')

% serial para paralelo
a_p = reshape(a, [8 length(a)/8]);

%% TX OFDM
sk = ifft(a_p, numSC);

figure()
hold on
stem(real(sk))
stem(imag(sk))
title('IFFT')
ylabel('s_k')
legend('real', 'imaginary')

txbb = sk; % TODO
pt = @(t) sqrt(T/numSC) * sin(pi*numSC*t/T) ./ (pi*t);
% st = conv (txbb, pt, 'same');
%SRF = st.*exp(j*2*pi*f*t)

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

rk = sk; % TODO
% y_p = sign(real(fft( rk, numSC )));
y_p = fft( rk, numSC );

size(y_p)

% paralelo para serial
y = reshape(y_p, [1 size(y_p, 1)*size(y_p, 2)]);

% slicer
z = (y >= 0) - (y < 0);

figure()
stem(z)
title('RX')
xlabel('samples')
ylabel('z')

%% RX BER
err = sum(double(a(1:length(z))) ~= z);
if (err >0 )
    error([num2str(err) ' errors'])
else
    disp('no errors')
end

