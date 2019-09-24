function [s, t] = mod_pam4(a)
%MOD_PAM4 Modulador para formato "Pulse Amplitude Modulation" (PAM) com
% 4-niveis

plot_en = true;

% frequencia de amostragem
fs = 8000; % Hz
% fs = 1000;

% frequencia central da portadora
% fc = 1000;
fc = 400; % Hz
% fc = 1;

len_sym = length(a); % quantidade de simbolos

R = fc;  % taxa de transmissao de simbolos
T = 1/R; % periodo de simbolo

%% Verificacoes de consistencia de configuracao
if fc >= fs/2
	error('A frequencia da portadora deve ser menos da metade da frequencia de amostragem.');
end

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

%% Tipos de Filtros 

% A) cosseno levantado
raised_cosine = @(t,a) sinc(t/T) .* (cos(a*pi*t/T) ./ (1 - (2*a*t/T).^2));

% B) onda quadrada
% square_wave = @(t) %TODO: implementar

%% Sinal formatador de pulso
alpha = 0.5;
gt = raised_cosine(t, alpha);

%% Pulsos formatados
s_conv = conv(am, gt);
% t2 = linspace(2*tmin, 2*tmax - (1/fs), length(s_conv));

t2_0 = t0;
s = s_conv(t2_0:t2_0+length(t)-1); % TODO: remover tempo anterior a zero

if (plot_en)
    figure()
    hold on
    % vetor de espacamentos em tempo discreto
    k = linspace(T, T * len_sym, len_sym);
    stem(k, a)
%     plot(t2, s_conv)
    plot(t, s)
    xlim([0 tmax])
    title('PAM-4 TX - Shaped Pulses')
    xlabel('Time')
    ylabel('s(t)')
    hold off
end

end

