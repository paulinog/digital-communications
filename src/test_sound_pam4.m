%% PAM-4 Transmitter
% gpaulino
clc; clearvars;
close all;
disp('PAM-4 example')
%% Extra information
% fc = 440
%
% sample rate
% fs = 8192

%% User parameters
% transmitter
len_sym = 80; % mapped symbols

symbols_set = [-3, -1, 1, 3];

alpha = 0; % raised cossine alpha

% time_interval = 10; % seconds
fs = 8192; % sampling rate (Hz)

fc = 1024; % carrier frequency
% R = fc; % symbol transmission rate
% T = 1/R; % symbol period

T = 1/fc;
R = 1/((1+alpha)*T);

% channel
% noise amplitude (Vpp ratio of the received signal)
% N0 = 0; 
N0 = 0.1;
% N0 = 0.5;

% general
plot_en = true;
plot_en_all = false;

%% Start a stopwatch timer
tic;

%% Time vector
time_interval = (len_sym +1) * (1/R);
tmin = -time_interval;
tmax = time_interval;

t = linspace(tmin, tmax, (tmax-tmin)*fs);
t0 = fs * (tmax-tmin)/2 + 1; % initial time

%% TX Source + Mapper
a = randsrc(1, len_sym, symbols_set); % random symbols source
m = linspace(t0 + (1/R)*fs, t0 + (1/R)*fs*len_sym, len_sym); % spacing vector continuous in time

am = zeros(1, length(t)); % spaced symbols
am(m) = a;

if (plot_en)
    % spacing vector discrete in time
    k = linspace((1/R), (1/R) * len_sym, len_sym);
    if (plot_en_all)
        figure()
        stem(k, a) 
        title('TX Random Symbols Source (4-levels)')
        xlabel('Samples')
        ylabel('Symbols')
        xlim([0 tmax])
    end
end
%% TX Raised cosine filter
g = @(t,a) sinc(t/T) .* (cos(a*pi*t/T) ./ (1 - (2*a*t/T).^2));

if (plot_en && plot_en_all)
    alphas = [0, 0.25, 0.5, 1]; % example figure
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

if (plot_en && plot_en_all)
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
    xlim([0 tmax])
    title('TX Shaped Pulses (PAM-4)')
    xlabel('Time')
    legend('a(k)', 's(t)')
    hold off
    movegui('northwest')
end
%% TX Up Converter

s_up = s .* sin(2*pi*fc .* t2);

if (plot_en)
    figure()
    hold on
    plot(t2,s)
    plot(t2,s_up)
    xlim([0 tmax])
    title('TX Up Converter (PAM-4)')
    xlabel('Time')
    legend('s(k)', 's\_up(t)')
    hold off
    movegui('northwest')
end

%% Channel
% N = N0*(2*max(s)); % noise amplitude
% %n = N*gaussmf(t2, [0.5 5]);
% n = N*rand(1, length(t2));
% 
% if (plot_en && plot_en_all)
%     figure()
%     plot(t2, n)
%     xlim([tmin tmax])
%     xlabel('Time')
%     ylabel('Amplitude')
%     title('White Noise')
% end
% %% RX Received Signals
% r = s+n;
% 
% if (plot_en && plot_en_all)
%     figure()
%     plot(t2,r)
%     xlim([min(t) max(t2)])
%     ylabel('r(t)')
%     xlabel('Time')
%     title('RX Received Signals')
% end

%% Sound

sound(s_up)
