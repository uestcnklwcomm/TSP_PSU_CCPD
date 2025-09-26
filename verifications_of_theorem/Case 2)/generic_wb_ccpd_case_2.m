%% generic case: proposed coupled CPD (single receiver)

clc; clear; close all;

Fs_nyquist = 2e9;         
Ts = 1 / Fs_nyquist; 
multicoset_N = 13;
Ns = 5e5;
Ns_all = Ns*multicoset_N;
t = (0:Ns_all-1)*Ts;  
undetermined = 1;  
rng(40, 'twister');
%% signal setup
R = 4;               % Number of transmissions
Bw = 5e7;
fc = [-5e8, 1e8, 5e8, -3e8, 4e8, 7e8]; %% carrier frequencies

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
Exp_linear = 0:1:(Mantenna-1);

H_gt = zeros(Mantenna,R);
for rr = 1:R
    multi_path = 10 + randperm(20,1);
    channel_gain = randn(1,multi_path) + 1j * randn(1,multi_path);
    channel_gain = channel_gain/(sum(channel_gain.*conj(channel_gain)));
    for pp = 1:multi_path
        theta_pp = -pi*sin(pi*rand);
        Exp_pp = exp(1j * Exp_linear * theta_pp);
        H_gt(:,rr) = H_gt(:,rr) + channel_gain(pp) * Exp_pp';
    end
end
y_nyquist = H_gt * signals; % Nyquist received sequence

%% multi-coset sampling 

sampled_positions = 0:multicoset_N:Ns_all-multicoset_N+1; %% multi-coset sampling
lag_channels = [1,2];

y_multicoset = cell(1,numel(lag_channels));
for ll = 1:numel(lag_channels)

    y_multicoset{ll} = y_nyquist(:, sampled_positions + lag_channels(ll));

end



freq_subNyquist = linspace(-Fs_nyquist/(2*multicoset_N),Fs_nyquist/(2*multicoset_N),num_fft);

%% cross-correlation

tensor_G0 = xcross_spectra(y_multicoset{1}, y_multicoset{1}, num_fft); %%cross-correlation tensor direct path


%% ccpd 
options_ft = 0;
Factor_ft = alg_ccpd_m3_fib(tensor_G0, R, options_ft);

xauto_spectra = Factor_ft{end};


%% estimating H (not necessary)
    
H_est = zeros(Mantenna,R);
m3_G0 = tens2mat(tensor_G0,[],3);

%% flipping negative PSDs, in case there's possibly a scaling
for rr = 1:R
    [absmax, pos] = max(abs(real(xauto_spectra(:,rr))));
    if absmax/real(xauto_spectra(pos,rr)) < 0
        xauto_spectra(:,rr) = -xauto_spectra(:,rr);
    end
end

krHHconj = m3_G0/(xauto_spectra.');

for rr = 1:R
    rank1r = krHHconj(:,rr);
    [H_est(:,rr),res] = rank_1_fac_conjugate(reshape(rank1r,[Mantenna,Mantenna]));
    % [Hestn(:,rr),res] = rank_1_phase_correction(reshape(rank1r,[Mantenna,Mantenna]));
end

Herr = cpderr(H_gt, H_est);


cross_spectra{1} = fftshift(xauto_spectra.',2);

%% cross-correlation tensor Gq
Ncomb = nchoosek(1:numel(lag_channels),2);
num_comb = size(Ncomb,1);


for cc = 1:num_comb

    receive_1 = y_multicoset{Ncomb(cc,1)};
    receive_2 = y_multicoset{Ncomb(cc,2)};
    
    tensor_Gq = xcross_spectra(receive_1, receive_2,num_fft ...
        );
    
    %% calibrating phase rotation 
    delta_t = lag_channels(Ncomb(cc,2)) - lag_channels(Ncomb(cc,1)); 
    time_shift = delta_t * Ts; 
    phase_correction = exp(1i * 2 * pi * freq_subNyquist * time_shift);
    m3_Gq = tens2mat(tensor_Gq,[],3);
    cross_spectra_cc = fftshift(pinv(krHHconj) * m3_Gq,2);
    cross_spectra{cc+1} = cross_spectra_cc.*phase_correction;

end


%% cross-correlation cpd
Ccol = zeros(num_comb+1,1);
Ccol(1) = 1;
%% creating phase lag
for cc = 1: num_comb
    lag_1 = lag_channels(Ncomb(cc,1));
    lag_2 = lag_channels(Ncomb(cc,2));
    Ccol(cc+1) = exp(1i*2*pi*(lag_1 - lag_2)/multicoset_N);
end
% 
Vandermonde_C = zeros(numel(Ccol),multicoset_N);

for mr = 1:multicoset_N
    Vandermonde_C(:,mr) = Ccol.^(mr-(multicoset_N+1)/2);
end

%% Establishing equations
alias_pxx_observation = cell(1,R);

for rr = 1:R
    
    alias_pxx_observation_r = zeros(1+num_comb,num_fft);

    for cc = 1:num_comb+1
        alias_pxx_observation_r(cc,:) = cross_spectra{cc}(rr,:);
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
        A_sel = Vandermonde_C(:,support_ff');
        spectra_est = A_sel \ alias_pxx_observation_r(:,ff);
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



