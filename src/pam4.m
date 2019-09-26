%% PAM-4 Transmitter
% gpaulino
clc; clearvars; close all;
%% Extra information
% fc = 440
%
% sample rate
% fs = 8192

%% User parameters
% transmitter
len_sym = 4; % mapped

symbols_set = [-3, -1, 1, 3];

alpha = 0.5; % raised cossine alpha

% time_interval = 10; % seconds
fs = 1000; % sampling rate

fc = 1; % carrier frequency
R = fc; % symbol transmission rate
T = 1/R; % symbol period

% channel
N0 = 0.1; % noise amplitude (Vpp ratio of the received signal)

% general
plot_en = true; % enable plot

%% Start a stopwatch timer
tic;

%% Time vector
time_interval = (len_sym +1) * T;
tmin = -time_interval;
tmax = time_interval;

t = linspace(tmin, tmax, (tmax-tmin)*fs);
t0 = fs * (tmax-tmin)/2 + 1; % initial time

%% TX Source + Mapper
a = randsrc(1, len_sym, symbols_set); % random symbols source
m = linspace(t0 + T*fs, t0 + T*fs*len_sym, len_sym); % spacing vector continuous in time

am = zeros(1, length(t)); % spaced symbols
am(m) = a;

if (plot_en)
    figure()
    % spacing vector discrete in time
    k = linspace(T, T * len_sym, len_sym);
    stem(k, a) 
    title('TX Random Symbols Source (4-levels)')
    xlabel('Samples')
    ylabel('Symbols')
    xlim([tmin tmax])
end
%% TX Raised cosine filter
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
%% TX Pulse Shape Filter
gt = g(t, alpha);

if (plot_en)
    figure()
    plot(t, gt)
    title(['TX Pulse Shape Filter - Raised cosine (\alpha=' num2str(alpha) ')'])
    xlabel('Time')
    ylabel('g(t)')
end
%% TX Shapped Pulses
s = conv(am, gt);
t2 = linspace(2*tmin, 2*tmax, 2*(tmax-tmin)*fs - 1);

if (plot_en)
    figure()
    hold on
    % stem(t,am)
    stem(k,a)
    plot(t2,s)
    xlim([min(t) max(t2)])
    title('TX Shaped Pulses (PAM-4)')
    xlabel('Time')
    ylabel('s(t)')
    hold off
end
%% Channel
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
%% RX Received Signals
r = s+n;

if (plot_en)
    figure()
    plot(t2,r)
    xlim([min(t) max(t2)])
    ylabel('r(t)')
    xlabel('Time')
    title('RX Received Signals')
end
%% RX Matched Filter
h = conj( g(-t, alpha) );

if (plot_en)
    figure()
    plot(t,h)
    title(['RX Matched Filter - Raised cosine (\alpha=' num2str(alpha) ')'])
    xlabel('Time')
    ylabel('h(t) = g*(-t)')
end
%% RX Normalize
y = conv(r, h);

t3 = linspace(3*tmin, 3*tmax, 3*(tmax-tmin)*fs - 2);
t0_3 = (((3*tmax)-(3*tmin))/2)*fs + 1; % sync point

norm_y = y;
% norm_y = y - abs(y(t0_3));
pos_y_ratio = abs(max(r)/max(y));
neg_y_ratio = abs(min(r)/min(y));

norm_y(norm_y > 0) = norm_y(norm_y > 0) .* pos_y_ratio;
norm_y(norm_y < 0) = norm_y(norm_y < 0) .* neg_y_ratio;

if (plot_en)
    figure()
    hold on
    stem(k,a)
    plot(t3, norm_y)
    title('RX Normalized Signal')
    xlabel('Time')
    legend('a(k)','norm\_y(t)')
    xlim([min(t) max(t2)])
    hold off
end
%% TX/RX Signals

if (plot_en)
    xmin = 0;
    xmax = 10;

    ymin = min([min(s) min(r) min(norm_y)]); % signal + noise
    ymax = max([max(s) max(r) max(norm_y)]);
    
    figure()
    subplot(411)
    stem(k,a)
    ylabel('a(k)')
    grid on
    xlim([xmin xmax])
    ylim([ymin ymax])
    title('1: Source Symbols')
    subplot(412)
    plot(t2,s)
    ylabel('s(t)')
    xlim([xmin xmax])
    ylim([ymin ymax])
    grid on
    title('2: Transmitted Signal')
    subplot(413)
    plot(t2, r)
    ylabel('r(t)')
    xlim([xmin xmax])
    ylim([ymin ymax])
    grid on
    title('3: Received Signal')
    subplot(414)
    plot(t3, norm_y)
    ylabel('norm\_y(t)')
    xlim([xmin xmax])
    ylim([ymin ymax])
    grid on
    title('4: Filtered Signal')
    xlabel('Time')
end
%% RX Quantization
m3 = linspace(1, len_sym, len_sym);
k3 = linspace(t0_3 + T*fs, 2*length(t)-fs+1, len_sym);  
fm = norm_y(k3); % sampled levels

if (plot_en)
    figure()
    hold on
    stem(m3, fm)
    plot(t3, norm_y)
    legend('f(k) = norm\_y(t=kT)','norm\_y(t)')
    xlim([0 max(t)])
    xlabel('Time')
    title('RX Quantization')
    hold off
end
%% RX Slicer
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
    stem(k, a)
    stem(m3, out)
    xlim([0 tmax])
    legend('a(k)', 'â(k)')
    title('RX Slicer')
    hold off
end
%%
disp(['input:  ' num2str(a)])
disp(['output: ' num2str(out)])
err = sum(ne(a,out));
if (err == 0)
    disp(['errors: ' num2str(err) ' out of total ' num2str(len_sym) ' symbols'])
else
    error(['errors: ' num2str(err) ' out of total ' num2str(len_sym) ' symbols'])
end
%% end section
toc;