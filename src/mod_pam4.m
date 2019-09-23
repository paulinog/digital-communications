function s = mod_pam4(a)
%MOD_PAM4 Modulador para formato "Pulse Amplitude Modulation" (PAM) com
% 4-niveis

%% User parameters
% transmitter
len_src_sym = 10; % mapped
T = 1;        % period in seconds
symbols_set = [-3, -1, 1, 3];

alpha = 0; % raised cossine alpha

% time_interval = 10; % seconds
samples_per_second = 1000; % samples

% channel
N0 = 0.1; % noise amplitude (Vpp ratio of the received signal)

% receiver
A_norm = 1.25; % normalization gain

% general
plot_en = true; % enable plot

%% Time vector
time_interval = len_src_sym +1;

tmin = -time_interval;
tmax = time_interval;
t = linspace(tmin, tmax, (tmax-tmin)*samples_per_second);
t0 = ((tmax-tmin)/2)*samples_per_second + 1; % initial time

%% TX Source + Mapper
% a = randsrc(1, len_src_sym, symbols_set); % random symbols source
m = linspace(1, len_src_sym, len_src_sym);  % spacing vector discrete in time
k = linspace(t0 + T*samples_per_second, length(t)-samples_per_second+1, len_src_sym);  % spacing vector continuous in time

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

end

