%% PAM-4 Transmitter
clc; clear all; close all;
%% User parameters
symbols = 9; % mapped
T = 1;        % period in seconds
symbols_set = [-3, -1, 1, 3];

alpha = 0; % raised cossine alpha

time_interval = 10; % seconds
samples_per_second = 1000; % samples

%% Time vector
tmin = -time_interval;
% tmin = 0;
tmax = time_interval;
t = linspace(tmin, tmax, (tmax-tmin)*samples_per_second);
if (tmax == -tmin)
%     t0 = ((tmax-tmin)/2)*samples_per_second; % initial time
    t0 = ((tmax-tmin)/2)*samples_per_second + 1; % initial time
elseif (tmin == 0)
    t0 = 1;
else
    % raise exception
    t0 = 1;
end

%% TX Source + Mapper
a = randsrc(1, symbols, symbols_set); % random symbols source
m = linspace(1, symbols, symbols);  % spacing vector discrete in time
k = linspace(t0 + 1*samples_per_second, length(t)-samples_per_second+1, symbols);  % spacing vector continuous in time

am = zeros(1, length(t)); % spaced symbols
am(k) = a;

figure()
% hold on %
stem(m, a)
% stem(t,am) %
title('4-PAM TX - Random Source')
xlabel('Samples')
ylabel('Symbols')
xlim([tmin tmax])
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

figure()
hold on
% stem(t,am)
stem(m,a)
plot(t2,s)
hold off
xlim([min(t) max(t2)])
%%  PAM-4 Receiver
N = 0.1*(2*max(s)); % noise amplitude
%n = N*gaussmf(t2, [0.5 5]);
n = N*rand(1, length(t2));
figure()
plot(t2, n)
%% 
r = s+n;
figure()
plot(t2,r)
xlim([min(t) max(t2)])
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
norm_y = y.*max(symbols_set)./max(y);

figure()
hold on
stem(m,a)
plot(t3, norm_y)
hold off
xlim([min(t) max(t2)])
%%
xmin = 0;
xmax = 10;

ymin = min([min(s) min(r) min(norm_y)]); % signal + noise
ymax = max([max(s) max(r) max(norm_y)]);

figure()
subplot(411)
stem(m,a)
grid on
xlim([xmin xmax])
ylim([ymin ymax])
title('1: Source Symbols')
subplot(412)
plot(t2,s)
xlim([xmin xmax])
ylim([ymin ymax])
grid on
title('2: Transmitted Signal')
subplot(413)
plot(t2, r)
xlim([xmin xmax])
ylim([ymin ymax])
grid on
title('3: Received Signal')
subplot(414)
plot(t3, norm_y)
xlim([xmin xmax])
ylim([ymin ymax])
grid on
title('4: Filtered Signal')

%% slicer
for R = 1:length(t3)
    if (mod(R,100) == 0)
        
    end
end
