function phase_offset = phase_compensation(x, fc, fs, timestep, k)
%PHASE COMPENSATION
sync_vec = double((mls(k, 1) > 0.5) - (mls(k, 1) <= 0.5));
len_r = length(x);
t = 0:timestep:(len_r-1)*timestep;
[num, den] = butter(10, fc*2*timestep, 'low');

%% course loop
center_phi = 0;
max_peak = 0;
for phi_course = 0 : pi/8 : 2*pi
    RRF = x .*cos(2*pi*fc*t + phi_course);
    LP = filtfilt(num, den, RRF) * 2;
    y = downsample(LP, fs);
    self_corr = xcorr(y, sync_vec);
    abs_pks = findpeaks(real(self_corr), 'NPeaks', 2,'SortStr','descend');
    if max(abs_pks) > max_peak
        max_peak = max(abs_pks);
        center_phi = phi_course;
    end
end

%% fine loop
points = 10;
phi_offset = center_phi;
for phi_fine = center_phi-pi/4 : pi/(2*points) : center_phi+pi/4
    RRF = x .*cos(2*pi*fc*t + phi_fine);
    LP = filtfilt(num, den, RRF) * 2;
    y = downsample(LP, fs);
    self_corr = xcorr(y, sync_vec);
    abs_pks = findpeaks(real(self_corr), 'NPeaks', 2,'SortStr','descend');
    if max(abs_pks) > max_peak
        max_peak = max(abs_pks);
        phi_offset = phi_fine;
    end
end

phase_offset = phi_offset;
end

