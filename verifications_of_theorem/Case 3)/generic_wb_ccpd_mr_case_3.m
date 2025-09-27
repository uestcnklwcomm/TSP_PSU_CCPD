%% generic case: proposed coupled CPD (Multi-receiver)
clc; clear; close all;

Fs_nyquist = 2e9;         
Ts = 1 / Fs_nyquist; 
multicoset_N = 21;
Ns = 5e5;
Ns_all = Ns*multicoset_N;
t = (0:Ns_all-1)*Ts;  
undetermined = 0;
rng(40, 'twister');
%% signal setup
R = 4;               % Number of transmissions
Bw = 1.6e8;            %% baseband bandwidth
fc = [-5e8, 1e8, 6e8, -3e8, 4e8, 7e8]; %% carrier frequency overlapped 



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

Dreceivers = 2; % two synchronized receivers

H_gt = cell(1,Dreceivers);
y_nyquist = cell(1, Dreceivers);
for dd = 1:Dreceivers
   
    H_gt_d = zeros(Mantenna,R);

    for rr = 1:R
        multi_path = 10 + randperm(10,1);
        channel_gain = randn(1,multi_path) + 1j * randn(1,multi_path);
        channel_gain = channel_gain/(sum(channel_gain.*conj(channel_gain)));
        for pp = 1:multi_path
            theta_pp = -pi*sin(pi*rand);
            Exp_pp = exp(1j * Exp_linear * theta_pp);
            H_gt_d(:,rr) = H_gt_d(:,rr) + channel_gain(pp) * Exp_pp';
        end
    end
    y_nyquist{dd} = H_gt_d * signals; % Nyquist received sequence
    H_gt{dd} = H_gt_d;
end




%% multi-coset sampling 
sampled_positions = 0:multicoset_N:Ns_all-multicoset_N+1; %% multi-coset sampling
lag_channels{1} = [1, 2, 4, 13];
lag_channels{2} = [1, 8, 13];

y_multicoset = cell(1,Dreceivers);

for dd = 1:Dreceivers

    for ll = 1:numel(lag_channels{dd})
    
        y_multicoset{dd}{ll} = y_nyquist{dd}(:, sampled_positions + lag_channels{dd}(ll));
    
    end

end


freq_subNyquist = linspace(-Fs_nyquist/(2*multicoset_N),Fs_nyquist/(2*multicoset_N),num_fft);

%% denoising       

noise_power_est = zeros(Dreceivers,Mantenna);

tensor_G0 = cell(1,Dreceivers);

for dd = 1:Dreceivers

    tensor_G0_d = xcross_spectra(y_multicoset{dd}{1}, y_multicoset{dd}{1}, num_fft); %%cross-correlation tensor direct path
    tensor_G0{dd} = tensor_G0_d;

end


%% ccpd 
options_ft = 1;
Factor_ft = alg_ccpd_mu_fib(tensor_G0, R, options_ft);
xauto_spectra = Factor_ft{end};

%% estimating kr(H, H*)
krHHconj = cell(1, Dreceivers);

for dd = 1:Dreceivers
    krHHconj{dd} = kr(Factor_ft{dd*2}, Factor_ft{2*dd-1});
end
krHH1conj_gt = kr(conj(H_gt{1}), H_gt{1});

factor1_err = cpderr(krHH1conj_gt, krHHconj{1});


%% flipping negative PSDs, in case there's possibly a scaling
for rr = 1:R
    [absmax, pos] = max(abs(real(xauto_spectra(:,rr))));
    if absmax/real(xauto_spectra(pos,rr)) < 0
        xauto_spectra(:,rr) = -xauto_spectra(:,rr);
        for dd = 1:Dreceivers
            krHHconj{dd}(:,rr) = -krHHconj{dd}(:,rr);
        end

    end
end


cross_spectra{1} = fftshift(xauto_spectra.',2);

num_Q = 2;
%% cross-correlation tensor Gq && establishing matrix C
Ccol(1) = 1;

for dd = 1:Dreceivers
    
    lag_comb = nchoosek(1:numel(lag_channels{dd}), 2);
    
    krHHconj_d = krHHconj{dd};

    for cc = 1:size(lag_comb,1)

        idx1 = lag_comb(cc, 1);
        idx2 = lag_comb(cc, 2);
        
        Seq1 = y_multicoset{dd}{idx1};
        Seq2 = y_multicoset{dd}{idx2};

        tensor_Gqd = xcross_spectra(Seq1, Seq2, num_fft);

        %% calibrating phase rotation 
        delta_t = lag_channels{dd}(idx2) - lag_channels{dd}(idx1); 
        time_shift = delta_t * Ts; 
        phase_correction = exp(1i * 2 * pi * freq_subNyquist * time_shift);

        Ccol(num_Q) = exp(1i*2*pi*(-delta_t)/multicoset_N);
        
        m3_Gqd = tens2mat(tensor_Gqd,[],3);
        cross_spectra{num_Q} = fftshift(pinv(krHHconj_d) * m3_Gqd,2) .*phase_correction ;
        num_Q = num_Q + 1;
    end

end


Vandermonde_C = zeros(num_Q-1,multicoset_N);

for mr = 1:multicoset_N
    Vandermonde_C(:,mr) = (Ccol.').^(mr-(multicoset_N+1)/2);
end

%% Establishing equations
alias_pxx_observation = cell(1,R);

for rr = 1:R
    
    alias_pxx_observation_r = zeros(num_Q-1,num_fft);

    for cc = 1:num_Q - 1
        alias_pxx_observation_r(cc,:) = cross_spectra{cc}(rr,:);
    end
    
    alias_pxx_observation{rr} = alias_pxx_observation_r;

end

%% PSD Unmixing
reconstructed_pxx = cell(1,R);
k_sparse = ceil(max(Bw)/(Fs_nyquist/multicoset_N));
power_spectra_est = zeros(multicoset_N * num_fft, R);

brute_force = 0; % 1: brute-force; 0: omp

for rr = 1:R
    
    reconstructed_pxx_r = zeros(num_fft,multicoset_N); %% size F*L, each row represents recovered PSD at L intervals

    alias_pxx_observation_r = alias_pxx_observation{rr};

    for ff = 1:num_fft
        if brute_force
            support_ff = bruteforce(alias_pxx_observation_r(:,ff), Vandermonde_C, k_sparse);
            C_sel = Vandermonde_C(:,support_ff');
            spectra_est = C_sel \ alias_pxx_observation_r(:,ff);
            reconstructed_pxx_r(ff,support_ff) = spectra_est.';
        else
            
            [spectra_est,~] = my_omp(Vandermonde_C, alias_pxx_observation_r(:,ff), k_sparse);
            reconstructed_pxx_r(ff,:) = spectra_est.';
        
        end

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



