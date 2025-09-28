clc;clear; close all;

load('Nov05-08.mat'); 
multicoset_N = 11;
num_fft = ceil(1024/multicoset_N); 
Mantennas = size(data,2);


Fs_nyquist = 3e5;
Ts = 1/Fs_nyquist;
sampled_position = 0:multicoset_N:size(data,1) - multicoset_N;
lag_channels = [1,2,4,8];
R = 2;

y_multicoset = cell(1,numel(lag_channels));
data_filtered = data.';

for ll = 1:numel(lag_channels)

    y_multicoset{ll} = data_filtered(:, sampled_position + lag_channels(ll));

end

freq_subNyquist = linspace(-Fs_nyquist/(2*multicoset_N),Fs_nyquist/(2*multicoset_N),num_fft);


%% denoising
coupled_tensor_G0 = xcross_spectra_nd(y_multicoset{1}, y_multicoset{1}, num_fft);
num_ccpd = length(coupled_tensor_G0);

%% ccpd 
options_ft = 1;
Factor_ft = alg_ccpd_mu_fib(coupled_tensor_G0, R, options_ft);

xauto_spectra = Factor_ft{end};
%% estimating H1H2
krHHconj = [];
for cc = 1:num_ccpd 
    
    hfac1 = Factor_ft{2*(cc-1)+1};
    hfac2 = Factor_ft{2*cc};
    hconj12 = kr(hfac2, hfac1);
    krHHconj = [krHHconj; hconj12];

end


%% flipping negative PSDs, in case there's possibly a scaling
for rr = 1:R
    [absmax, pos] = max(abs(real(xauto_spectra(:,rr))));
    if absmax/real(xauto_spectra(pos,rr)) < 0
        xauto_spectra(:,rr) = -xauto_spectra(:,rr);
        krHHconj(:,rr) = -krHHconj(:,rr);
    end
end

cross_spectra{1} = fftshift(xauto_spectra.',2);

%% cross-correlation tensor Gq
Ncomb = nchoosek(1:numel(lag_channels),2);
num_comb = size(Ncomb,1);


for cc = 1:num_comb

    receive_1 = y_multicoset{Ncomb(cc,1)};
    receive_2 = y_multicoset{Ncomb(cc,2)};
    
    coupled_tensor_Gq = xcross_spectra_nd(receive_1, receive_2,num_fft ...
        );
    
    m3_Gq = [];
    
    for ccc = 1:num_ccpd
            
        m3_Gq = [m3_Gq; tens2mat(coupled_tensor_Gq{ccc},[],3)];
    
    end


    %% calibrating phase rotation 
    delta_t = lag_channels(Ncomb(cc,2)) - lag_channels(Ncomb(cc,1)); 
    time_shift = delta_t * Ts; 
    phase_correction = exp(1i * 2 * pi * freq_subNyquist * time_shift);

    cross_spectra_cc = fftshift(pinv(krHHconj) * m3_Gq,2);
    cross_spectra{cc+1} = cross_spectra_cc.*phase_correction;

end

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
    if mod(multicoset_N,2)
        Vandermonde_C(:,mr) = Ccol.^(mr-(multicoset_N+1)/2);
    else
        Vandermonde_C(:,mr) = Ccol.^(mr-(multicoset_N+2)/2);
    end
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
k_sparse = 3;

% [~, freq_plot_upsamp] = pwelch(signal_wb_total, hamming(num_fft*multicoset_L), [], [], Fs_nyquist, 'centered');

for rr = 1:R
    
    reconstructed_pxx_r = zeros(num_fft,multicoset_N); %% size F*L, each row represents recovered PSD at L intervals
    alias_pxx_observation_r = alias_pxx_observation{rr};
    % support_rr = support_set{rr};
    for ff = 1:num_fft
    
        support_ff = bruteforce(alias_pxx_observation_r(:,ff), Vandermonde_C, k_sparse);
        C_sel = Vandermonde_C(:,support_ff);
        spectra_est = C_sel \ alias_pxx_observation_r(:,ff);
        reconstructed_pxx_r(ff,support_ff) = spectra_est.';

    end

    PSD_est_normalized(rr,:) = abs(reconstructed_pxx_r(:).')/max(abs( ...
        reconstructed_pxx_r(:).'));
end




support_1 = [round(205/1024 * num_fft * multicoset_N):round(410/1024 * num_fft * multicoset_N), ...
    round(615/1024 * num_fft * multicoset_N): round(819/1024 * num_fft * multicoset_N)];
support_2 = [round(411/1024 * num_fft *multicoset_N):round(614/1024*num_fft*multicoset_N)];
SMR_proposed = smr_measure(PSD_est_normalized(1,:), support_1, PSD_est_normalized(2,:), support_2);

mixed_spectra = xcross_spectra(data_filtered, data_filtered, 1024); 
mixed_spectra_3 = fftshift(abs(squeeze(mixed_spectra(3,3,:))) ...
    /(max(abs(squeeze(mixed_spectra(3,3,:))))));

figure(1); clf;

tfig = tiledlayout(2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

custom_colors = [0 0 1;     % blue
                 1 0 0;     % red
                 1 1 0;     % yellow
                 0.5 0 0.5];% purple
nexttile;

plot(1: 1024, mixed_spectra_3, ...
        'Color', custom_colors(1,:), ...
        'LineWidth', 1); 
set(gca, 'FontName', 'Times New Roman', 'FontSize', 14);
xlabel('Frequency bins', 'FontName', 'Times New Roman', 'FontSize', 14);
xlim([1 1024]);
ylabel('Normalized PSDs', 'FontName', 'Times New Roman', 'FontSize', 14);

grid on;
set(gca, 'GridLineStyle', '--', 'LineWidth', 1);


title('Mixed PSDs', ...
      'Interpreter', 'latex', ...
      'FontName', 'Times New Roman', ...
      'FontSize', 14);
legend off;


nexttile;
for rr = 1:R  
    plot(1: num_fft *multicoset_N, PSD_est_normalized(rr,:), ...
        'Color', custom_colors(rr,:), ...
        'LineWidth', 1); 
    hold on;
end


set(gca, 'FontName', 'Times New Roman', 'FontSize', 14);
xlabel('Frequency bins', 'FontName', 'Times New Roman', 'FontSize', 14);
xlim([1 num_fft*multicoset_N]);
ylabel('Normalized PSDs', 'FontName', 'Times New Roman', 'FontSize', 14);

grid on;
set(gca, 'GridLineStyle', '--', 'LineWidth', 1);


title(['Separated PSDs, SMR = ' , sprintf('%0.4g', SMR_proposed)], ...
      'Interpreter', 'latex', ...
      'FontName', 'Times New Roman', ...
      'FontSize', 14);

legend off;


