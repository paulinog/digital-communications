%% PAM-4 Transmitter
clc; clear all; close all;
%% User parameters
symbols = 10; % mapped
T = 1;        % period in seconds
symbols_set = [-3, -1, 1, 3];

alpha = 0; % raised cossine alpha

time_interval = 10; % seconds
samples_per_second = 1000; % samples

%% Time vector
% tmin = -time_interval;
tmin = 0;
tmax = time_interval;
t = linspace(tmin, tmax, (tmax-tmin)*samples_per_second);
%% TX Source + Mapper
a = randsrc(1, symbols, symbols_set); % random symbols source
m = linspace(0, symbols-1, symbols);  % spacing vector discrete in time
% k = linspace(0, length(t)/2-samples_per_second, symbols)+1;  % spacing vector continuous in time
k = linspace(0, length(t)-samples_per_second, symbols)+1;  % spacing vector continuous in time

am = zeros(1, length(t)); % spaced symbols
% am = zeros(1, length(t)/2); % spaced symbols
am(k) = a;

figure()
% hold on
stem(m, a)
% stem(t,am)
title('4-PAM TX - Random Source')
xlabel('Samples')
ylabel('Symbols')
%% Raised cosine filter
g = @(t,a) sinc(t/T) .* (cos(a*pi*t/T) ./ (1 - (2*a*t/T).^2));

alphas = [0, 0.5, 0.8, 1]; % example figure
figure()
subplot(221)
plot(t, g(t, alphas(1)))
title(['Raised cosine \alpha=' num2str(alphas(1))])
ylabel('g(t)')
subplot(222)
plot(t, g(t, alphas(2)))
title(['Raised cosine \alpha=' num2str(alphas(2))])
subplot(223)
plot(t, g(t, alphas(3)))
title(['Raised cosine \alpha=' num2str(alphas(3))])
xlabel('Time')
subplot(224)
plot(t, g(t, alphas(4)))
title(['Raised cosine \alpha=' num2str(alphas(4))])

%% Pulse Shape Filter
gt = g(t, alpha);
figure()
plot(t, gt)
title(['Raised cosine (\alpha=' num2str(alpha) ')'])
xlabel('Time')
ylabel('g(t)')
%% Shapped Pulses
s = conv(am, gt);
t2 = linspace(2*tmin, 2*tmax, 2*(tmax-tmin)*samples_per_second - 1);
% t2 = linspace(2*tmin, 2*tmax, 2*(tmax-tmin/2)*samples_per_second - 1);

figure()
hold on
plot(t2,s)
stem(m,a)
% stem(t,am)
hold off
%%  PAM-4 Receiver
N = 0.5*max(s);
%n = N*gaussmf(t2, [0.5 5]);
n = N*rand(1, length(t2));
figure(4)
plot(t2, n)
%% 
r = s+n;
figure()
plot(t2,r)
%%
h = conj( g(-t, alpha) );
figure()
plot(t,h)
title(['Matched Filter - Raised cosine (\alpha=' num2str(alpha) ')'])
xlabel('Time')
ylabel('g*(-t)')
%% 
t3 = linspace(3*tmin, 3*tmax, 3*(tmax-tmin)*samples_per_second - 2);
% t3 = linspace(3*tmin, 3*tmax, 4*(tmax-tmin/2)*samples_per_second - 1);
y = conv(r, h);
figure()
plot(t3, y)
%%
figure()
subplot(411)
stem(m,a)
title('Source Symbols')
subplot(412)
plot(t2,s)
title('Transmitted Signal')
subplot(413)
plot(t2, r)
title('Received Signal')
subplot(414)
plot(t3, y)
title('Filtered Signal')

%% slicer
for R = 1:length(t3)
    if (mod(R,100) == 0)
        
    end
end
