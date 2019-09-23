function s = mod_pam4(a)
%MOD_PAM4 Modulador para formato "Pulse Amplitude Modulation" (PAM) com
% 4-niveis

alpha=0.5;
plot_en = true;

% samples_per_second = 8000;
% fc = 400;

samples_per_second = 1000;
fc = 1;


len_src_sym = length(a);
T = 1/fc;
%% Time vector
time_interval = len_src_sym +1;

tmin = -time_interval;
tmax = time_interval;
t = linspace(tmin, tmax, (tmax-tmin)*samples_per_second);
t0 = ((tmax-tmin)/2)*samples_per_second + 1;

%%
m = linspace(1, len_src_sym, len_src_sym);  % spacing vector discrete in time
k = linspace(t0 + T*samples_per_second, length(t)-samples_per_second+1, len_src_sym);  % spacing vector continuous in time

am = zeros(1, length(t)); % spaced symbols
am(k) = a;

%% Raised cosine filter
g = @(t,a) sinc(t/T) .* (cos(a*pi*t/T) ./ (1 - (2*a*t/T).^2));

%% Pulse Shape Filter
gt = g(t, alpha);

%% Shapped Pulses
s = conv(am, gt);
t2 = linspace(2*tmin, 2*tmax, 2*(tmax-tmin)*samples_per_second - 1);

if (plot_en)
    figure()
    hold on
    % stem(t,am)
    stem(m, a)
    plot(t2,s)
    xlim([min(t) max(t2)])
    title(['4-PAM TX - Shaped Pulses'])
    xlabel('Time')
    ylabel('s(t)')
    hold off
end

end

