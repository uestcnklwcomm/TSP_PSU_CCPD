%% generic case: SPA
clc; clear; close all;


Fs_nyquist = 2e6; %Nyquist rate 
Ts = 1 / Fs_nyquist;       
Ns = 5e5;                 
t = (0:Ns-1)*Ts;            
noiseless = 1;
snr = 5;
rng(40,'twister');
%% setup
R = 4;               % Number of transmissions
Bw = 3e5;
fc = [-7e5, -1e5, 3e5, 7e5];  % central frequency of each transmission

signals = zeros(R,Ns);

% baseband white gaussian noise
white_Gaussian = (randn(R, Ns) + 1i * randn(R, Ns)) / sqrt(2); 
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

%% generating CSI

Mantenna = 4; %needs over-determined systems
Exp_linear = 0:1:(Mantenna-1);

H_gt = zeros(Mantenna,R);
for rr = 1:R
    multi_path = 10 + randperm(10,1);
    channel_gain = randn(1,multi_path) + 1j * randn(1,multi_path);
    channel_gain = channel_gain/(sum(channel_gain.*conj(channel_gain)));
    for pp = 1:multi_path
        theta_pp = -pi*sin(pi*rand);
        Exp_pp = exp(1j * Exp_linear * theta_pp);
        H_gt(:,rr) = H_gt(:,rr) + channel_gain(pp) * Exp_pp';
    end
end
y_nyquist = H_gt * signals; % Nyquist received sequence

%% adding noise

if ~noiseless
    snr_linear = 10^(snr/10);
    for mm = 1:Mantenna
    
        energy = mean(y_nyquist(mm,:) .* conj(y_nyquist(mm,:)));
        energy_noise = energy/snr_linear;
        wg_noise = sqrt(energy_noise)* (randn(1, Ns ...
            ) + 1j * randn(1, Ns));
        hann_window = hann(Ns).';
        window_power = mean(abs(hann_window).^2);
        scale_factor = sqrt(Ns) / sqrt(window_power * Ns);  
        y_nyquist(mm,:) = y_nyquist(mm,: ...
            ) + scale_factor * hann_window .* wg_noise;
        
    end
end

%% NMF-SPA

NMF_mat = zeros(Mantenna,num_fft);

%% cross-correlation & denoising
for mm = 1:Mantenna

    NMF_mat(mm,:) = fftshift(abs(powerspectra_gt(y_nyquist(mm,:), num_fft)));

    if ~noiseless
        
        NMF_mat(mm,:) = NMF_mat(mm,:) - min(NMF_mat(mm,:));

    end

end

sel_ind = SPA(NMF_mat,R);
tilde_H = NMF_mat(:,sel_ind);
power_spectra_est = (pinv(tilde_H) * NMF_mat).';

% %% fine-tunning via NMF-HALS
% 
% max_iter = 1000;
% [~, power_spectra_est] = NMF_HALS(NMF_mat, R, max_iter, tilde_H, power_spectra_est);


freq = linspace(-Fs_nyquist/2, Fs_nyquist/2, num_fft);

figure(1); clf;

tfig = tiledlayout(2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

custom_colors = [0 0 1;     % blue
                 1 0 0;     % red
                 1 1 0;     % yellow
                 0.5 0 0.5];% purple

nexttile;
for rr = 1:R  
    Pxx = powerspectra_gt(signals(rr,:), num_fft);
    plot(freq, fftshift(abs(Pxx) ./ max(abs(Pxx))), ...
        'Color', custom_colors(rr,:), ...
        'LineWidth', 1); 
    hold on;
    P_gt(:,rr) = fftshift(abs(Pxx)).';
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

NMSE = cpderr(P_gt, power_spectra_est);

nexttile;
for rr = 1:R
    power_spectra_est_r = power_spectra_est(:, rr);
    plot(freq, (abs(power_spectra_est_r) ./ max(abs(power_spectra_est_r))), ...
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


