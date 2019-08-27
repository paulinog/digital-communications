%% PAM-4 Transmitter
% gpaulino
clc; clear all; close all;
%% User parameters
% transmitter
symbols = 9; % mapped
T = 1;        % period in seconds
symbols_set = [-3, -1, 1, 3];

alpha = 0; % raised cossine alpha

time_interval = 10; % seconds
samples_per_second = 1000; % samples

% channel
N0 = 0.1; % noise amplitude (Vpp ratio of the received signal)

% receiver
A_norm = 1.25; % normalization gain

% general
plot_en = true; % enable plot

%% user parameter conscistency self-check
assert(symbols*T < time_interval, 'number of symbols must fit in the time interval')


%% Time vector
tmin = -time_interval;
tmax = time_interval;
t = linspace(tmin, tmax, (tmax-tmin)*samples_per_second);
t0 = ((tmax-tmin)/2)*samples_per_second + 1; % initial time

%% TX Source + Mapper
a = randsrc(1, symbols, symbols_set); % random symbols source
m = linspace(1, symbols, symbols);  % spacing vector discrete in time
k = linspace(t0 + T*samples_per_second, length(t)-samples_per_second+1, symbols);  % spacing vector continuous in time

am = zeros(1, length(t)); % spaced symbols
am(k) = a;

if (plot_en)
    figure()
    % hold on %
    stem(m, a)
    % stem(t,am) %
    title('4-PAM TX - Random Source')
    xlabel('Samples')
    ylabel('Symbols')
    xlim([tmin tmax])
end
%% Raised cosine filter
g = @(t,a) sinc(t/T) .* (cos(a*pi*t/T) ./ (1 - (2*a*t/T).^2));

if (plot_en)
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
end
%% Pulse Shape Filter
gt = g(t, alpha);

if (plot_en)
    figure()
    plot(t, gt)
    title(['Raised cosine (\alpha=' num2str(alpha) ')'])
    xlabel('Time')
    ylabel('g(t)')
end
%% Shapped Pulses
s = conv(am, gt);
t2 = linspace(2*tmin, 2*tmax, 2*(tmax-tmin)*samples_per_second - 1);

if (plot_en)
    figure()
    hold on
    % stem(t,am)
    stem(m,a)
    plot(t2,s)
    xlim([min(t) max(t2)])
    title('4-PAM RX - Shaped Pulses')
    xlabel('Time')
    ylabel('s(t)')
    hold off
end
%%  PAM-4 Receiver
N = N0*(2*max(s)); % noise amplitude
%n = N*gaussmf(t2, [0.5 5]);
n = N*rand(1, length(t2));

if (plot_en)
    figure()
    plot(t2, n)
    xlim([tmin tmax])
    xlabel('Time')
    ylabel('Amplitude')
    title('White Noise')
end
%% 
r = s+n;

if (plot_en)
    figure()
    plot(t2,r)
    xlim([min(t) max(t2)])
    xlabel('Time')
    title('Received Signal')
end
%%
h = conj( g(-t, alpha) );

if (plot_en)
    figure()
    plot(t,h)
    title(['Matched Filter - Raised cosine (\alpha=' num2str(alpha) ')'])
    xlabel('Time')
    ylabel('g*(-t)')
end
%% normalize
y = conv(r, h);
norm_y = y.*A_norm*max(symbols_set)./max(abs(y));

t3 = linspace(3*tmin, 3*tmax, 3*(tmax-tmin)*samples_per_second - 2);
t0_3 = (((3*tmax)-(3*tmin))/2)*samples_per_second + 1; % sync point

if (plot_en)
    figure()
    hold on
    stem(m,a)
    plot(t3, norm_y)
    xlim([min(t) max(t2)])
    hold off
end
%%
xmin = 0;
xmax = 10;

ymin = min([min(s) min(r) min(norm_y)]); % signal + noise
ymax = max([max(s) max(r) max(norm_y)]);

if (plot_en)
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
end
%% quantization
m3 = linspace(1, symbols, symbols);
k3 = linspace(t0_3 + T*samples_per_second, 2*length(t)-samples_per_second+1, symbols);  
fm = norm_y(k3); % sampled levels

if (plot_en)
    figure()
    hold on
    stem(m3, fm)
    plot(t3, norm_y)
    % plot(norm_y)
    xlim([0 max(t)])
    hold off
end
%% slicer
out = zeros(1, length(fm));
i=1;

for R = fm
    if R > 2
        out(i) = 3;
    elseif R > 0
        out(i) = 1;
    elseif R > -2
        out(i) = -1;
    else
        out(i) = -3;
    end
    i = i + 1;
end

if (plot_en)
    figure()
    hold on %
    stem(m, a)
    stem(m3, out)
    hold off
end
%%
disp(['input:  ' num2str(a)])
disp(['output: ' num2str(out)])
err = sum(ne(a,out));
disp(['errors: ' num2str(err) ' out of total ' num2str(symbols) ' symbols'])