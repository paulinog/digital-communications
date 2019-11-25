function [phase_offset, sample_shift] = phase_compensation(x, fc, fs, timestep, k)
%PHASE COMPENSATION
sync_vec = double((mls(k, 1) > 0.5) - (mls(k, 1) <= 0.5));
% sync_bits = 2^k-1;
len_x = length(x);
t = 0:timestep:(len_x-1)*timestep;
[num, den] = butter(10, fc*2*timestep, 'low');

%% finding peaks 
% optimize the time of execution

% RRF_I = x(t_shift:end) .*cos(2*pi*fc*t(t_shift:end));
% RRF_Q = x(t_shift:end) .* -sin(2*pi*fc*t(t_shift:end));
% LP_I = filtfilt(num, den, RRF_I) * 2;
% LP_Q = filtfilt(num, den, RRF_Q) * 2;
% LP = LP_I + 1j*LP_Q;
% y = downsample(LP, fs);
% self_corr = xcorr(y, sync_vec);
% [abs_pks, abs_loc] = findpeaks(abs(self_corr), 'NPeaks', 2,'SortStr','descend');
% start_y = min(abs_loc) + sync_bits + 1 - length(y)
t_opt = t;

%% course shift
% initial sample shift
max_shift = 1;
points_course = 10;
max_peak = 0;
for shift_course = ceil(1 : fs/points_course : fs)
    RRF_I = x(shift_course:end) .*cos(2*pi*fc*t(shift_course:end));
    RRF_Q = x(shift_course:end) .* -sin(2*pi*fc*t(shift_course:end));
    LP_I = filtfilt(num, den, RRF_I) * 2;
    LP_Q = filtfilt(num, den, RRF_Q) * 2;
    LP = LP_I + 1j*LP_Q;
    y = downsample(LP, fs);
    self_corr = xcorr(y, sync_vec);
    abs_pks = findpeaks(abs(self_corr), 'NPeaks', 2,'SortStr','descend');
    if max(abs_pks) > max_peak
        max_peak = max(abs_pks);
        max_shift = shift_course;
    end
end
%% fine shift
points_fine = 10;
range_fine = fs/points_course;
step_fine = range_fine/(points_fine/2);
end_fine = max_shift + range_fine;
if max_shift > 1
    init_fine = max_shift - range_fine;
else
    init_fine = 1;
end
for shift_fine = ceil(init_fine:step_fine:end_fine)
    RRF_I = x(shift_fine:end) .*cos(2*pi*fc*t(shift_fine:end));
    RRF_Q = x(shift_fine:end) .* -sin(2*pi*fc*t(shift_fine:end));
    LP_I = filtfilt(num, den, RRF_I) * 2;
    LP_Q = filtfilt(num, den, RRF_Q) * 2;
    LP = LP_I + 1j*LP_Q;
    y = downsample(LP, fs);
    self_corr = xcorr(y, sync_vec);
    abs_pks = findpeaks(abs(self_corr), 'NPeaks', 2,'SortStr','descend');
    if max(abs_pks) > max_peak
        max_peak = max(abs_pks);
        max_shift = shift_fine;
    else
        break;
    end
end
disp('sample shift = ')
disp(max_shift)
t_sh = t_opt(max_shift:end);
% plot self corr, after time shift
% RRF_I = x(max_shift:end) .*cos(2*pi*fc*t(max_shift:end));
% RRF_Q = x(max_shift:end) .* -sin(2*pi*fc*t(max_shift:end));
% LP_I = filtfilt(num, den, RRF_I) * 2;
% LP_Q = filtfilt(num, den, RRF_Q) * 2;
% LP = LP_I + 1j*LP_Q;
% y = downsample(LP, fs);
% self_corr = xcorr(y, sync_vec);
% figure()
% plot(abs(self_corr));
% title('time shift')
% ylabel('R')
% xlabel('sample')

%% course loop
center_phi = 0;
phi_step = pi/8;
max_peak = 0;
for phi_course = -pi : phi_step : pi
    RRF = x(max_shift:end) .*cos(2*pi*fc*t_sh + phi_course);
    LP = filtfilt(num, den, RRF) * 2;
    y = downsample(LP, fs);
    self_corr = xcorr(y, sync_vec);
    pks_course = findpeaks(real(self_corr(length(y):end)), 'NPeaks', 2,'SortStr','descend');
    if max(pks_course) > max_peak
        max_peak = max(pks_course);
        center_phi = phi_course;
    end
end
% plot self corr, before fine loop
% RRF_I = x .*cos(2*pi*fc*t + center_phi);
% RRF_Q = x .* -sin(2*pi*fc*t + center_phi);
% LP_I = filtfilt(num, den, RRF_I) * 2;
% LP_Q = filtfilt(num, den, RRF_Q) * 2;
% LP = LP_I + 1j*LP_Q;
% y = downsample(LP, fs);
% self_corr = xcorr(y, sync_vec);
% figure()
% hold on
% plot(abs(self_corr));
% plot(real(self_corr));
% plot(imag(self_corr));
% title('center \phi')
% ylabel('R')
% xlabel('sample')
% legend('|R|','Re(R)','Im(R)')
% hold off

%% fine loop
phi_offset = center_phi;
points = 10;
for phi_fine = center_phi-phi_step : 2*phi_step/points : center_phi+phi_step
    RRF = x(max_shift:end) .*cos(2*pi*fc*t_sh + phi_fine);
    LP = filtfilt(num, den, RRF) * 2;
    y = downsample(LP, fs);
    self_corr = xcorr(y, sync_vec);
    pks_fine = findpeaks(real(self_corr(length(y):end)), 'NPeaks', 2,'SortStr','descend');
    if max(pks_fine) > max_peak
        max_peak = max(pks_fine);
        phi_offset = phi_fine;
    end
end
disp('phase offset (degrees) = ')
disp(phi_offset*360/(2*pi))

sample_shift = max_shift;
phase_offset = phi_offset;
end