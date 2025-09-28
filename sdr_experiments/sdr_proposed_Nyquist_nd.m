clc;clear; close all;

load('Nov05-08.mat'); 
num_fft = 1024; 
Mantennas = size(data,2);


Fs_nyquist = 3e5;
Ts = 1/Fs_nyquist;
R = 2;

data_filtered = data.';
y_nyquist = data_filtered;

coupled_tensor_G = xcross_spectra_nd(y_nyquist, y_nyquist, num_fft);
options_ft = 1;
Factor_ft = alg_ccpd_mu_fib(coupled_tensor_G , R, options_ft);

PSD_est_normalized = fftshift(abs(Factor_ft{end}.'),2);
PSD_est_normalized = eye(R)/diag(max(PSD_est_normalized.') ...
    )*PSD_est_normalized;

% for rr = 1:R
%     PSD_est_normalized(rr,:) = 
% end

support_1 = [205:410, 615:819];
support_2 = 411:614;
SMR_proposed = smr_measure(PSD_est_normalized(1,:), support_2, ...
    PSD_est_normalized(2,:), support_1);


mixed_spectra = xcross_spectra(data_filtered, data_filtered, num_fft); 
mixed_spectra_3 = fftshift(abs(squeeze(mixed_spectra(3,3,:))) ...
    /(max(abs(squeeze(mixed_spectra(3,3,:))))));

figure(1); clf;

tfig = tiledlayout(2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

custom_colors = [0 0 1;     % blue
                 1 0 0;     % red
                 1 1 0;     % yellow
                 0.5 0 0.5];% purple
nexttile;

plot(1: num_fft, mixed_spectra_3, ...
        'Color', custom_colors(1,:), ...
        'LineWidth', 1); 
set(gca, 'FontName', 'Times New Roman', 'FontSize', 14);
xlabel('Frequency bins', 'FontName', 'Times New Roman', 'FontSize', 14);
xlim([1 num_fft]);
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
    plot(1: num_fft, PSD_est_normalized(rr,:), ...
        'Color', custom_colors(rr,:), ...
        'LineWidth', 1); 
    hold on;
end


set(gca, 'FontName', 'Times New Roman', 'FontSize', 14);
xlabel('Frequency bins', 'FontName', 'Times New Roman', 'FontSize', 14);
xlim([1 num_fft]);
ylabel('Normalized PSDs', 'FontName', 'Times New Roman', 'FontSize', 14);

grid on;
set(gca, 'GridLineStyle', '--', 'LineWidth', 1);


title(['Separated PSDs, SMR = ' , sprintf('%0.4g', SMR_proposed)], ...
      'Interpreter', 'latex', ...
      'FontName', 'Times New Roman', ...
      'FontSize', 14);

legend off;


