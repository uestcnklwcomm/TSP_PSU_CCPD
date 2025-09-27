clc; clear; close all;


Fs_nyquist = 2e6;        
Ts = 1 / Fs_nyquist;       
Ns_all = 1e6;                 
t = (0:Ns_all-1)*Ts;         
snr = 5;
noiseless = 0;
undetermined = 0;
%% signal setup
rng(40, 'twister');
R = 4;               % Number of transmissions
Bw = 3e5;
fc = [-7e5, -1e5, 3e5, 7e5];  % central frequency of each transmission

signals = zeros(R,Ns_all);

% baseband white gaussian noise
white_Gaussian = (randn(R, Ns_all) + 1i * randn(R, Ns_all)) / sqrt(2); 
[Unitary, ~] = qr(randn(R) + 1j * randn(R));
Stats = Unitary * white_Gaussian;

for rr = 1:R
    fc_r = fc(rr);

    % low-pass filter
    f_cutoff = Bw / 2;
    filter_order = 300;
    lp = fir1(filter_order, f_cutoff/(Fs_nyquist/2), 'low', hamming(filter_order+1));
    x_bb = filter(lp, 1, Stats(rr,:));  % baseband WSS

    x_rf = x_bb .* exp(1j*2*pi*fc_r*t);

    signals(rr, :) = x_rf;
end

num_fft = 1024;
freq_Nyquist = linspace(-Fs_nyquist/2, Fs_nyquist/2, num_fft);


%% generating CSI
Mantenna = undetermined*3 + ~undetermined*4;
DOA_R = randperm(180,R); %DOAs of sources
DOA_R_d = -pi * sin(DOA_R/180*pi);% half wavelength distance
phi_row = exp(1i*DOA_R_d);

for mm = 1:Mantenna
    Phi(mm,:) = phi_row.^(mm-1);
end

a_vec = randn(1,R) + 1j * randn(1,R);
A = diag(a_vec);

y_nyquist = Phi * A * signals; % Nyquist received sequence


%% adding noise 

if ~noiseless
    snr_linear = 10^(snr/10);
    for mm = 1:Mantenna
    
        energy = mean(y_nyquist(mm,:) .* conj(y_nyquist(mm,:)));
        energy_noise = energy/snr_linear;
        wg_noise = sqrt(energy_noise)* (randn(1, Ns_all ...
            ) + 1j * randn(1, Ns_all));
        hann_window = hann(Ns_all).';
        window_power = mean(abs(hann_window).^2);
        scale_factor = sqrt(Ns_all) / sqrt(window_power * Ns_all);
        y_nyquist(mm,:) = y_nyquist(mm,: ...
            ) + scale_factor * hann_window .* wg_noise;
        
    end
end

%% denoising
denoising = 1;
tensor_G = xcross_spectra_LOS(y_nyquist, y_nyquist, R, num_fft, denoising);
tensor_G = permute(tensor_G, [1,3,2]);


%% unmixing
options_ft = 1;
Factor_ft = alg_ccpd_m3_fib(tensor_G, R, options_ft);
PSD_est = abs(Factor_ft{end});


figure(1); clf;

tfig = tiledlayout(2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

custom_colors = [0 0 1;     % blue
                 1 0 0;     % red
                 1 1 0;     % yellow
                 0.5 0 0.5];% purple

nexttile;
for rr = 1:R 
    Pxx = powerspectra_gt(signals(rr,:), num_fft);
    plot(freq_Nyquist, fftshift(abs(Pxx) ./ max(abs(Pxx))), ...
        'Color', custom_colors(rr,:), ...
        'LineWidth', 1); 
    hold on;
    P_gt(:,rr) = abs(Pxx).';
end


set(gca, 'FontName', 'Times New Roman', 'FontSize', 14);
xlabel('Frequency (Hz)', 'FontName', 'Times New Roman', 'FontSize', 14);
ylabel('Normalized PSDs', 'FontName', 'Times New Roman', 'FontSize', 14);

grid on;
set(gca, 'GridLineStyle', '--', 'LineWidth', 1);


title('Ground-truth, $R=4$, $M=3$', ...
      'Interpreter', 'latex', ...
      'FontName', 'Times New Roman', ...
      'FontSize', 14);

legend off;

NMSE = cpderr(P_gt, PSD_est);

nexttile;
for rr = 1:R
    PSD_estr = PSD_est(:, rr);
    plot(freq_Nyquist, fftshift(abs(PSD_estr) ./ max(abs(PSD_estr))), ...
        'Color', custom_colors(rr,:), ...
        'LineWidth', 1);
    hold on;
end
grid on;
set(gca, 'GridLineStyle', '--', 'LineWidth', 1);
xlabel('Frequency (Hz)', 'FontName', 'Times New Roman', 'FontSize', 14);
ylabel('Normalized PSDs', 'FontName', 'Times New Roman', 'FontSize', 14);
title(['Estimated PSDs, NMSE = ', sprintf('%0.4g', NMSE)], 'Interpreter', 'latex', ...
      'FontName', 'Times New Roman', ...
      'FontSize', 14);
set(gca, 'FontName', 'Times New Roman', 'FontSize', 14);
legend off;

