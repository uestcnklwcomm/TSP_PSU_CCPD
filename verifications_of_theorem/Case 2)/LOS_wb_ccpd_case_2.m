%% generic case: proposed coupled CPD (single receiver)

clc; clear; close all;

Fs_nyquist = 2e9;         
Ts = 1 / Fs_nyquist;   
multicoset_N = 13;
% Ns = 5e5;
Ns_all = 1e6;
t = (0:Ns_all-1)*Ts;
undetermined = 1;
rng(40, 'twister');
%% signal setup
R = 4;               % Number of transmissions
Bw = 5e7;
fc = [-5e8, 1e8, 5e8, -3e8, 4e8, 7e8]; %% carrier frequency overlapped 
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

num_fft = 64;

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



%% multi-coset sampling 
sampled_positions = 0:multicoset_N:Ns_all-multicoset_N+1; % multi-coset sampling
lag_channels = [1,7];

y_multicoset = cell(1,numel(lag_channels));
for ll = 1:numel(lag_channels)

    y_multicoset{ll} = y_nyquist(:, sampled_positions + lag_channels(ll));

end

freq_subNyquist = linspace(-Fs_nyquist/(2*multicoset_N),Fs_nyquist/(2*multicoset_N),num_fft);


Qcomb = nchoosek(1:numel(lag_channels),2);
Q = size(Qcomb,1);

%% cross-correlation of direct path      
tensor_G{1} = xcross_spectra_LOS(y_multicoset{1}, y_multicoset{1}, R, num_fft, 0);

%% delayed path
for qq = 1:Q

    Seq1 = y_multicoset{Qcomb(qq,1)};
    Seq2 = y_multicoset{Qcomb(qq,2)};
    tensor_G{qq+1} = xcross_spectra_LOS(Seq1, Seq2, R, num_fft, 0);

end


%% ccpd 
options_ft = 1;
Factor_ft = alg_ccpd_mu_fib(tensor_G, R, options_ft);

%% estimating aliased power spectra by unifying the ambiguities
krHPhi1 = Factor_ft{1};
cross_spectra{1} = fftshift(Factor_ft{2}.',2);

for qq = 1:Q

    krHPhin = Factor_ft{2*qq+1};
    [err,~,scaling] = cpderr(krHPhi1, krHPhin); 
    cross_spectra_q = Factor_ft{2*(qq+1)}/scaling{1};
    delta_t = lag_channels(Qcomb(qq,2)) - lag_channels(Qcomb(qq,1)); 
    time_shift = delta_t * Ts; 
    phase_correction = exp(1i * 2 * pi * freq_subNyquist * time_shift);
    cross_spectra{qq+1} = fftshift(cross_spectra_q.',2).*phase_correction;

end


Ccol = zeros(Q+1,1);
Ccol(1) = 1;
%% creating phase lag
for qq = 1: Q

    lag_1 = lag_channels(Qcomb(qq,1));
    lag_2 = lag_channels(Qcomb(qq,2));
    Ccol(qq+1) = exp(1i*2*pi*(lag_1 - lag_2)/multicoset_N);

end
% 
Vandermonde_C = zeros(numel(Ccol),multicoset_N);

for mr = 1:multicoset_N
    Vandermonde_C(:,mr) = Ccol.^(mr-(multicoset_N+1)/2);
end

%% Establishing equations
alias_pxx_observation = cell(1,R);

for rr = 1:R
    
    alias_pxx_observation_r = zeros(Q+1,num_fft);

    for qq = 1:Q+1
        alias_pxx_observation_r(qq,:) = cross_spectra{qq}(rr,:);
    end
    
    alias_pxx_observation{rr} = alias_pxx_observation_r;

end

%% PSD Unmixing
reconstructed_pxx = cell(1,R);
k_sparse = ceil(max(Bw)/(Fs_nyquist/multicoset_N));
power_spectra_est = zeros(multicoset_N * num_fft, R);

for rr = 1:R
    
    reconstructed_pxx_r = zeros(num_fft,multicoset_N); %% size F*L, each row represents recovered PSD at L intervals

    alias_pxx_observation_r = alias_pxx_observation{rr};

    for ff = 1:num_fft
    
        support_ff = bruteforce(alias_pxx_observation_r(:,ff), Vandermonde_C,k_sparse);
        C_sel = Vandermonde_C(:,support_ff');
        spectra_est = C_sel \ alias_pxx_observation_r(:,ff);
        reconstructed_pxx_r(ff,support_ff) = spectra_est.';

    end

    power_spectra_est(:,rr) = abs(reconstructed_pxx_r(:));

end

%% ground-truth
power_spectra_gt = zeros(num_fft*multicoset_N, R);
for rr = 1:R  
    Pxx = powerspectra_gt(signals(rr,:), num_fft*multicoset_N);
    power_spectra_gt(:,rr) = fftshift(abs(Pxx)).';
end

nmse = cpderr(power_spectra_gt, power_spectra_est);

       

freq_Nyquist = linspace(-Fs_nyquist/2, Fs_nyquist/2, num_fft*multicoset_N);

figure(1); clf;

tfig = tiledlayout(2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

custom_colors = [0 0 1;     % blue
                 1 0 0;     % red
                 1 1 0;     % yellow
                 0.5 0 0.5];% purple

nexttile;
for rr = 1:R 
    power_spectra_gt_r = power_spectra_gt(:,rr).';
    plot(freq_Nyquist, power_spectra_gt_r/max(power_spectra_gt_r), ...
        'Color', custom_colors(rr,:), ...
        'LineWidth', 1); 
    hold on;
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


nexttile;
for rr = 1:R
    power_spectra_estr = power_spectra_est(:, rr).';
    plot(freq_Nyquist,  power_spectra_estr/ max(power_spectra_estr), ...
        'Color', custom_colors(rr,:), ...
        'LineWidth', 1);
    hold on;
end
grid on;
set(gca, 'GridLineStyle', '--', 'LineWidth', 1);
xlabel('Frequency (Hz)', 'FontName', 'Times New Roman', 'FontSize', 14);
ylabel('Normalized PSDs', 'FontName', 'Times New Roman', 'FontSize', 14);
title(['Estimated PSDs, NMSE = ', sprintf('%0.4g', nmse)], 'Interpreter', 'latex', ...
      'FontName', 'Times New Roman', ...
      'FontSize', 14);
set(gca, 'FontName', 'Times New Roman', 'FontSize', 14);
legend off;



