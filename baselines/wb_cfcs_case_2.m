%% Wide-band: CF-CS (This kind of approach usually need over-determined CPD system to accurately estimate the carrier frequecies)
clc; clear; close all;

Fs_nyquist = 2e9;         
Ts = 1 / Fs_nyquist;       
Ns = 5e5;                 
            
R = 2;
Bw = 4.5e7;

fc = [4e8, 1e8, 5e8, 6e8]; %carrier frequencies

delay_bw = (1/Bw);
delay_fc = (1/max(abs(fc)));
Mantenna = 4;
sub_Nyquist_rate = 13;
Ns_all = Ns * sub_Nyquist_rate;

t = (0:Ns_all-1)*Ts;

snr = 5;
noiseless = 0;


rng(40, 'twister')
delays = rand(1, Mantenna-1) * min(delay_bw, delay_fc);
%% signal generation

white_Gaussian = (randn(R, Ns_all) + 1i * randn(R, Ns_all)) / sqrt(2); 
[Unitary, ~] = qr(randn(R) + 1j * randn(R));
Stats = Unitary * white_Gaussian;

signals  = zeros(R,Ns_all);

for rr = 1:R
    fc_r = fc(rr);

    % low-pass filter
    f_cutoff = Bw / 2;
    filter_order = 300;
    lp = fir1(filter_order, f_cutoff/(Fs_nyquist/2), 'low', hamming(filter_order+1));
    x_bb{rr} = filter(lp, 1, Stats(rr,:));  % baseband WSS

    x_rf = x_bb{rr} .* exp(1j*2*pi*fc_r*t);

    signals(rr, :) = x_rf;
end

num_fft = 128; %% number of fft per coset channel

%% generating CSI
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
y_nyquist = H_gt * signals;

if ~noiseless
    snr_linear = 10^(snr/10);
    
    for mm = 1:Mantenna
    
        energy = mean(y_nyquist(mm,:) .* conj(y_nyquist(mm,:)));
        energy_noise = energy/snr_linear;
        wg_noise = sqrt(energy_noise)* (randn(1, Ns_all) + 1j * randn(1, Ns_all));
        hann_window = hann(Ns_all).';
        window_power = mean(abs(hann_window).^2);
        scale_factor = sqrt(Ns_all) / sqrt(window_power * Ns_all);
    
        y_nyquist(mm,:) = y_nyquist(mm,:) + scale_factor * hann_window .* wg_noise;
        
    end
end


y_delay = zeros(Mantenna, Ns_all);
y_delay(1, :) = y_nyquist(1, :);

for mm = 2: Mantenna
    
    for rr = 1:R  
        fc_r = fc(rr);
        signals_delay(rr,:) = x_bb{rr} .* exp(1j*2*pi*fc_r*(t + delays(mm - 1)));
    end

    y_delay(mm,:) = H_gt(mm, :) * signals_delay;

end

%% adding noise for delay path
if ~noiseless

    for mm = 2:Mantenna
    
        energy = mean(y_delay(mm,:) .* conj(y_delay(mm,:)));
        energy_noise = energy/snr_linear;
        wg_noise = sqrt(energy_noise)* (randn(1, Ns_all) + 1j * randn(1, Ns_all));
        hann_window = hann(Ns_all).';
        window_power = mean(abs(hann_window).^2);
        scale_factor = sqrt(Ns_all) / sqrt(window_power * Ns_all);
    
        y_delay(mm,:) = y_delay(mm,:) + scale_factor * hann_window .* wg_noise;
        
    end
end






sampled_positions = 0: sub_Nyquist_rate: Ns_all - sub_Nyquist_rate +1;
y_direct_channel = y_nyquist(:, sampled_positions + 1);
y_delay_channel = y_delay(:, sampled_positions + 1);




%% formulating CCPD
tensor_G_direct = zeros(Mantenna, Mantenna, 2*num_fft-1);
tensor_G_delay = zeros(Mantenna, Mantenna, 2*num_fft-1);


for mm = 1:Mantenna
    
    seq_direct_1 = y_direct_channel(mm,:);
    seq_delay_1 = y_delay_channel(mm,:);

    for mmm = 1:Mantenna

        seq_direct_2 = y_direct_channel(mmm,:);

        seq_delay_2 = y_delay_channel(mmm,:);

        tensor_G_direct(mm,mmm,:) =  (1/(2*num_fft-1)) * xcorr(seq_direct_1.',seq_direct_2.',num_fft-1).';

        tensor_G_delay(mm,mmm,:) = (1/(2*num_fft-1)) * xcorr(seq_delay_1.',seq_delay_2.',num_fft-1).';

    end
end

%% solving cpd via als
U0{1} = randn(Mantenna,R);
U0{2} = randn(Mantenna,R);
U0{3} = randn(2*num_fft - 1,R);

options.MaxIter = 500;
options.TolX = 1e-6;

%% direct path
[U_direct,~] = cpd_als(tensor_G_direct, U0, options);

%% re-initialize
U0{1} = randn(Mantenna,R);
U0{2} = randn(Mantenna,R);
U0{3} = randn(2*num_fft - 1,R);

%% delay path`
[U_delay,~] = cpd_als(tensor_G_delay, U0, options);


%% calibrating permutation mismatch
Auto_corr_direct = U_direct{end};
Auto_corr_delay = U_delay{end};

[Err,permute,~,~] = cpderr(Auto_corr_direct, Auto_corr_delay);

Hest_direct = U_direct{1};
Hest_delay = U_delay{1} * permute;

scaling = diag(Hest_direct(1,:)./Hest_delay(1,:));

Hest_delay_calibrate = Hest_delay * scaling;

%% estimating carrier frequency
fc_est = zeros(1,R);

for rr = 1:R
    omega_r_est = 0;
    for mm = 2: Mantenna

        hmr_direct = Hest_direct(mm,rr);
        hmr_delay = Hest_delay_calibrate(mm,rr);

        angle_direct = angle(hmr_direct);

        angle_delay = angle(hmr_delay);
        
%         residual = 
        omega_r_est = omega_r_est + mod(angle_delay - angle_direct, 2*pi)/(delays(mm-1));

    end
   
    fc_est(rr) = omega_r_est/(Mantenna-1)/(2*pi);
    
%     fc_est(rr) = omega_r_est/(Mantenna-1);
end


%% reconstructing power spectra
freq = linspace(-Fs_nyquist/2, Fs_nyquist/2, num_fft * sub_Nyquist_rate);
freq_start = freq(1: num_fft : num_fft*(sub_Nyquist_rate -1) + 1);
freq_end = freq(num_fft : num_fft : num_fft * sub_Nyquist_rate);


power_spectra_est = zeros(R, num_fft * sub_Nyquist_rate);


for rr = 1:R
    
    power_spectra_est_r = zeros(num_fft, sub_Nyquist_rate);

    psd_r_baseband = zeros(1, num_fft);

    for ff = 1:num_fft
                
        psd_r_baseband(ff) = (exp(-1i*2*pi*((ff-1)/num_fft)*[-(num_fft-1):(num_fft-1)])*Auto_corr_direct(:,rr))/num_fft;
    
    end

    interval_idx = find(freq_start <= fc_est(rr) & fc_est(rr) <= freq_end);

    power_spectra_est_r(:,interval_idx) = abs(fftshift(psd_r_baseband)).';

    power_spectra_est(rr,:) = (power_spectra_est_r(:)).';
        
end

% freq_upsamp = linspace(-Fs/)

figure(1); clf;

tfig = tiledlayout(2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

custom_colors = [0 0 1;     % blue
                 1 0 0;     % red
                 1 1 0;     % yellow
                 0.5 0 0.5];% purple

nexttile;
for rr = 1:R  
    Pxx = powerspectra_gt(signals(rr,:), num_fft*sub_Nyquist_rate);
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


title('Ground-truth, $R=2$, $M=4$', ...
      'Interpreter', 'latex', ...
      'FontName', 'Times New Roman', ...
      'FontSize', 14);

legend off;

NMSE = cpderr(P_gt, power_spectra_est.');

nexttile;
for rr = 1:R
    power_spectra_est_r = power_spectra_est(rr,:);
    plot(freq, power_spectra_est_r/max(power_spectra_est_r), ...
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

