function s = mod_pam4(a)
%MOD_PAM4 Modulador para formato "Pulse Amplitude Modulation" (PAM) com
% 4-niveis

alpha=0.5;
plot_en = true;

% frequencia de amostragem
fs = 8000; % Hz
% fs = 1000;

% frequencia central da portadora
% fc = 1000;
fc = 400; % Hz
% fc = 1;

len_sym = length(a);
T = 1/fc;

%% Vetor tempo
time_interval = (len_sym +1) * T;
tmin = -time_interval;
tmax = time_interval;

t = linspace(tmin, tmax, (tmax-tmin)*fs);
t0 = fs * (tmax-tmin)/2 + 1;

%% Sequencia espacada em tempo continuo
m = linspace(t0 + T*fs, t0 + T*fs*len_sym, len_sym);
am = zeros(1, length(t)); % simbolos espacados
am(m) = a;

%% Filtro cosseno levantado
g = @(t,a) sinc(t/T) .* (cos(a*pi*t/T) ./ (1 - (2*a*t/T).^2));

%% Sinal formatador de pulso
gt = g(t, alpha);

%% Pulsos formatados
s = conv(am, gt);
t2 = linspace(2*tmin, 2*tmax - (1/fs), length(s));

if (plot_en)
    figure()
    hold on
    % vetor de espacamentos em tempo discreto
    k = linspace(T, T * len_sym, len_sym);
    stem(k, a)
    plot(t2, s)
    xlim([0 tmax])
    title('PAM-4 TX - Shaped Pulses')
    xlabel('Time')
    ylabel('s(t)')
    hold off
end

end

