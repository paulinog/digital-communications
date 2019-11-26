%% Projeto Final - Parte 2
% Disciplina IE533A, 2s2019 (FEEC/UNICAMP)
% Guilherme Paulino, RA 117119
clc;
clearvars;
close all;

disp('Final Project IE533 - Part 2')
disp('Acoustic OFDM Modulator')

%% Parametros de usuario
N = 8; % N-point FFT
fc = 440; % frequencia da portadora, em Hz
fs = 8192; % frequencia de amostragem, em Hz
T = 1.3653; % periodo do simbolo, em segundos
% R = N/T;

% gravador de audio
audio_enable = true;
% true: *.wav file, false: loopback
audio_file_nloopback = false; 

% mostra todas as figuras
enable_plot = false;

%% Parametros do FEC
enable_FEC = true;
trellis = poly2trellis(3, [5 7]);
parity_ratio = 2;

%% Parametros do MLS
enable_MLS = true;
k = 5;
sync_bits = 2^k-1;

%% Fonte
% 1) Entrada de texto pelo usuario
% disp('Enter text:');
% input_str = input('','s');

% 2) Entrada de teste
input_str = 'abcdefghijklmnopqrstuvxzwyABCDEFGHIJKLMNOPQRSTUVXZWY1234567890';

%% Vetor binario
input_bin = str_source(input_str);

%% Mapeador
a = mapper(input_bin, [0 1]);
numBits = length(a);
disp(['Number of bits sent:' num2str(numBits)]);

if enable_plot
    figure()
    stem(a)
    title('TX')
    xlabel('samples')
end

%% TX Vetor tempo
timestep = T/(N*fs);

frame_size = numBits;
if enable_FEC
    frame_size = frame_size * parity_ratio;
end
if enable_MLS
    frame_size = frame_size + 2*sync_bits;
end
num_padding_t = zero_padding(frame_size, N);
if num_padding_t > 0
    frame_size = frame_size + num_padding_t;
end
tmax = frame_size*T/N;
t = 0:timestep:tmax-timestep;
t_pt = -tmax/2:timestep:tmax/2-timestep;

%% TX FEC
if enable_FEC
    a_enc = convenc(a, trellis);
else
    a_enc = a;
end

%% TX BPSK
a_mod = (a_enc >= 0.5) - (a_enc < 0.5);

%% TX Serial para paralelo
num_padding_mod = zero_padding(length(a_mod), N);
if num_padding_mod > 0
    a_mod = [a_mod zeros(1, num_padding_mod)];
end
an = reshape(a_mod, [N length(a_mod)/N]);

%% TX OFDM Mux
skn = ifft(an, N);

if enable_plot
    figure()
    subplot(211)
    stem(abs(skn)')
    title('IFFT')
    ylabel('|s_k|')
    subplot(212)
    stem(angle(skn)')
    title('IFFT')
    ylabel('\angle{s_k}')
    xlabel(['samples / ' num2str(N)])
end

%% TX Paralelo para serial
sk = reshape(skn, [1 size(skn, 1)*size(skn, 2)]);

%% TX MLS
if enable_MLS
    sync_vec = (1/2)*double((mls(k, 1) > 0.5) - (mls(k, 1) <= 0.5));
    sk_sync = [sync_vec sk sync_vec]; % concatenate
    
    num_padding_tx_mls = zero_padding(length(sk_sync), N);
    if num_padding_tx_mls > 0
        sk_sync = [sk_sync zeros(1, num_padding_tx_mls)];
    end
else
    sk_sync = sk;
end

%% TX Upsample
sk_up = upsample(sk_sync, fs, fs/2);
t_tx = 0:timestep:(length(sk_up)-1)*timestep;

if enable_plot
    figure()
    subplot(211)
    plot(t_tx,abs(sk_up))
    ylabel('|s_{up}[k]|')
    title('s[k] upsampled')
    subplot(212)
    plot(t_tx,angle(sk_up))
    xlabel('time')
    ylabel('\angle{s_{up}[k]}')
    xlabel('time')
end 

%% TX Formatador de pulso
% p = @(t) sqrt(T/N) * sin(pi*N*t/T) ./ (pi*t);
p = @(t) sqrt(N/T) * sinc(N*t/T);
pt = p(t_pt);

if enable_plot
    figure();
    plot(t_pt, pt);
    xlabel('time')
    ylabel('p(t)')
    title('Pulse shapper function')
end

%% TX OFDM Banda base
st = conv (sk_up, pt, 'same');

if enable_plot
    % amplitude e fase
    figure()
    subplot(211)
    plot(t_tx,abs(st))
    xlabel('time')
    ylabel('|s(t)|')
    subplot(212)
    plot(t_tx,angle(st))
    xlabel('time')
    ylabel('\angle{s(t)}')
    % sinal complexo
    figure()
    hold on
    plot(t_tx, real(st))
    plot(t_tx, imag(st))
    legend('real','imag')
    xlabel('time')
    ylabel('s(t)')
end

%% TX OFDM Banda passante
SRF = real(st).*cos(2*pi*fc*t_tx) - imag(st).*sin(2*pi*fc*t_tx);

if 1%enable_plot
    figure()
    hold on
    plot(t_tx, real(SRF))
    plot(t_tx, imag(SRF))
    legend('real','imag')
    xlabel('time')
    ylabel('SRF')
end

%% Audio Recorder
if audio_enable
    numBits = 16;
    numChannels = 1;
    Fs = round(1/timestep);
    
    disp('transmit ->')
    if audio_file_nloopback
        audio_rx = audiorecorder(Fs, numBits, numChannels);
        audiowrite('audio_file.wav', SRF/16, Fs)
        input('Press ENTER to start recording')
        record(audio_rx);
        input('Press ENTER to stop')
    else
        audio_rx = audiorecorder(Fs, numBits, numChannels);
        tic
        disp('Start recording audio')
        record(audio_rx);
        sound(SRF/16, Fs, numBits)
        pause(t(end) + 0.1); % delay
        toc
    end
    disp('receive <-')
    stop(audio_rx)
    disp('Processing...')
    
    r = 9*getaudiodata(audio_rx)';
else
    %% Channel
    % % noiseless
    % r = SRF;

    % % real valued noise
    % N0 = 0.8*max(st);
    % n = N0*rand(1, length(t));
    % r = SRF + n;

    % complex noise
    % snr = -40; %dB
    snr = -5; %dB
    % snr = 0;
    powerDB = 10*log10(var(SRF));
    noiseVar = 10.^(0.1*(powerDB-snr)); 
    r = awgn(SRF, snr);

    % TEST channel delays
    ch_delay1 = round(100*numBits*rand(1)); % random spacing
    ch_delay2 = round(100*numBits*rand(1));
    % r = [zeros(1, ch_delay1) r];
    r = [zeros(1, ch_delay1) r zeros(1, ch_delay2)];
end

%% RX Vetor tempo
len_r = length(r);
t_rx = 0:timestep:(len_r-1)*timestep;

if 1%enable_plot 
    figure()
    plot(t_rx, r)
    title('Received signal')
    ylabel('r(t) = s(t) + n(t)')
    xlabel('time')
end 

%% RX RF
if enable_MLS
    tic
    [phaseRF, sample_shift] = phase_compensation(r, fc, fs, timestep, k);
    disp('Time to compute phase compensation:')
    toc
    t_rx = t_rx(sample_shift:end);
else
    phaseRF = 0;
    sample_shift = 1;
end
RRF_I = r(sample_shift:end) .*cos(2*pi*fc*t_rx + phaseRF);
RRF_Q = r(sample_shift:end) .* -sin(2*pi*fc*t_rx + phaseRF);

if enable_plot
    figure()
    hold on
    plot(t_rx, RRF_I)
    plot(t_rx, RRF_Q)
    legend('real','imag')
    xlabel('time')
    ylabel('RRF')
end 

%% RX LP
[num, den] = butter(10, fc*2*timestep, 'low');
LP_I = filtfilt(num, den, RRF_I) * 2;
LP_Q = filtfilt(num, den, RRF_Q) * 2;

if enable_plot
    figure()
    hold on
    plot(t_rx, LP_I)
    plot(t_rx, LP_Q)
    legend('real','imag')
    xlabel('time')
    ylabel('LP')
end

%% RX Complex signal
rt = LP_I + 1j*LP_Q;

if enable_plot
    figure()
    subplot(211)
    plot(t_rx,abs(rt))
    xlabel('time')
    ylabel('|r(t)|')
    subplot(212)
    plot(t_rx,angle(rt))
    xlabel('time')
    ylabel('\angle{r(t)}')
end

%% RX OFDM downsample
% rk_up = rt(t = kT);
rk = downsample(rt, fs, fs/2);

%% RX Remover MLS
if enable_MLS
    sync_vec2 = (1/2)*double((mls(k, 1) > 0.5) - (mls(k, 1) <= 0.5));
    self_corr = xcorr(rk, sync_vec2);
    
    [pks, loc] = findpeaks(abs(self_corr), 'NPeaks', 2,'SortStr','descend');
    start_frame = min(loc) + sync_bits + 1 - length(rk);
    end_frame = max(loc) - length(rk);
    frame_size = end_frame - start_frame + 1;
    
    if 1%enable_plot
        figure()
        hold on
        plot(abs(self_corr));
        plot(real(self_corr));
        plot(imag(self_corr));
        title('Cross correlation')
        ylabel('R')
        xlabel('sample')
        legend('|R|','Re(R)','Im(R)')
        hold off
    end
else
    start_frame = 1;
    end_frame = length(rk);
    frame_size = end_frame - start_frame + 1;
end
rk_enc = rk(start_frame : end_frame);

% disp('R_start =')
% disp(start_frame)
% disp('R_end =')
% disp(end_frame)
% disp('frame size =')
% disp(frame_size)

%% RX Serial para paralelo
num_padding_rx_mls = zero_padding(length(rk_enc), N);
if num_padding_rx_mls > 0
    rk_enc = [rk_enc zeros(1, num_padding_rx_mls)];
end
rkn = reshape(rk_enc, [N length(rk_enc)/N]);

%% RX OFDM Demux
% yn = sign(real(fft( rkn, N )));
yn = fft( rkn, N );

%% RX Paralelo para serial
y = reshape(yn, [1 size(yn, 1)*size(yn, 2)]);

%% RX Slicer
z_sync = (y > 0);

if enable_plot
    figure()
    stem(z_sync)
    title('RX')
    xlabel('samples')
    ylabel('z')
end

%% RX FEC
if enable_FEC
    z = vitdec(z_sync, trellis, 20,  'trunc', 'hard');
else
    z = z_sync;
end

%% Demapeador
output_bin = demapper(z, [0 1]);

%% RX Destino
output_str = str_dest(output_bin);
% Mostra a saida de texto
disp('-----')
disp('Output text:');
fprintf('%s', output_str)

%% BER
disp(' ')
disp('-----')
% Calcula a taxa de erro de bit
if length(z) >= length(a)
    err = sum(a ~= z(1:length(a)));
    if (err > 0)
       if (err < 10)
           loc_err = find(a ~= z(1:length(a)));
           disp(['Position: ' num2str(loc_err)])
           disp(['Expected: ' num2str(a(loc_err))])
           disp(['Received: ' num2str(z(loc_err))])
       end
       
       error(['Total of errors: ' num2str(err) ' (' ...
               num2str(100*err/frame_size) '%)']);
    else
        disp('no errors');
    end
else
    error('wrong size of z')
end
