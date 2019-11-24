function phase_offset = phase_compensation(x, fc, fs, timestep, k)
%PHASE COMPENSATION
sync_vec = double((mls(k, 1) > 0.5) - (mls(k, 1) <= 0.5));
sync_bits = 2^k-1;
len_x = length(x);
t = 0:timestep:(len_x-1)*timestep;
[num, den] = butter(10, fc*2*timestep, 'low');

%% finding peaks
% RRF_I = x .*cos(2*pi*fc*t);
% RRF_Q = x .* -sin(2*pi*fc*t);
% LP_I = filtfilt(num, den, RRF_I) * 2;
% LP_Q = filtfilt(num, den, RRF_Q) * 2;
% LP = LP_I + 1j*LP_Q;
% y = downsample(LP, fs);
% self_corr = xcorr(c, sync_vec);
% [abs_pks, abs_loc] = findpeaks(abs(self_corr), 'NPeaks', 2,'SortStr','descend');
% start_y = min(abs_loc) + sync_bits + 1 - length(y)

t_opt = t;

%% course loop
phi_step = pi/8;
center_phi = 0;
max_peak = 0;
for phi_course = -pi : phi_step : pi
    RRF = x .*cos(2*pi*fc*t_opt + phi_course);
    LP = filtfilt(num, den, RRF) * 2;
    y = downsample(LP, fs);
    self_corr = xcorr(y, sync_vec);
    pks_course = findpeaks(real(self_corr(length(y):end)), 'NPeaks', 2,'SortStr','descend');
    if max(pks_course) > max_peak
        max_peak = sum(pks_course);
        center_phi = phi_course;
    end
end

RRF_I = x .*cos(2*pi*fc*t + center_phi);
RRF_Q = x .* -sin(2*pi*fc*t + center_phi);
LP_I = filtfilt(num, den, RRF_I) * 2;
LP_Q = filtfilt(num, den, RRF_Q) * 2;
LP = LP_I + 1j*LP_Q;
y = downsample(LP, fs);
self_corr = xcorr(y, sync_vec);
figure()
hold on
plot(abs(self_corr));
plot(real(self_corr));
plot(imag(self_corr));
title('center \phi')
ylabel('R')
xlabel('sample')
legend('|R|','Re(R)','Im(R)')
hold off

%% fine loop
points = 50;
phi_offset = center_phi;
for phi_fine = center_phi-phi_step : 2*phi_step/points : center_phi+phi_step
    RRF = x .*cos(2*pi*fc*t_opt + phi_fine);
    LP = filtfilt(num, den, RRF) * 2;
    y = downsample(LP, fs);
    self_corr = xcorr(y, sync_vec);
    pks_fine = findpeaks(real(self_corr(length(y):end)), 'NPeaks', 2,'SortStr','descend');
    if max(pks_fine) > max_peak
        max_peak = sum(pks_fine);
        phi_offset = phi_fine;
    end
end

phase_offset = phi_offset;
end

