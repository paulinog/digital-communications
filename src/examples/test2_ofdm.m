close all;
clc; clearvars; 
disp('OFDM example')
%% User parameters
numSymbol = 100;
fc0 = 440;
fc1 = 4*fc0;
fs = 4*fc1;

%% Time vector
timestep = 1/fs;
T = 1/fc0;

frame_size = numSymbol;
tmax = frame_size * T;

t = 0:timestep:tmax-timestep;

%% TX
msg = randsrc(1, numSymbol, [-1 1]);

figure()
stem(msg)
title('TX')
xlabel('samples')

%% TX OFDM
a = ifft( msg );
sk = a;
st = sk;
%% Channel 
rt = st;

figure()
hold on
stem(real(rt))
stem(imag(rt))
hold off
title('Complex Channel')
ylabel('r(n)')
xlabel('samples')
legend('real', 'imaginary')

% plot(t_rx, r)
% ylabel('r(t)')
% xlabel('time')

%% RX OFDM
rk = rt;
y = sign(real(fft( rk )));

figure()
stem(y)
title('RX')
xlabel('samples')
ylabel('y')

%% RX BER
err = sum(msg ~= y);
if (err >0 )
    error([num2str(err) ' errors'])
else
    disp('no errors')
end

