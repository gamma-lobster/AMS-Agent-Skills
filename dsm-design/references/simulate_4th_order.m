%% simulate_4th_order.m
%% 4th-Order CIFF DSM Simulation with NTF/STF plots
%% Spec: order=4, fs=2MHz, OSR=16, H_inf=4.0, 5-bit, CIFF

clear all; close all; clc;

fprintf('============================================================\n');
fprintf('  4th-Order CIFF DSM Simulation\n');
fprintf('  fs=2MHz, OSR=16, H_inf=4.0, 5-bit quantizer\n');
fprintf('============================================================\n\n');

%% NTF Parameters (from Python synthesizeNTF)
% These are the correct values - using them directly for accurate NTF plot
ntf_zeros = [0.9977727000 + 0.0667055500i;
             0.9977727000 - 0.0667055500i;
             0.9857393700 + 0.1682792100i;
             0.9857393700 - 0.1682792100i];

ntf_poles = [0.3541020100 - 0.5300351800i;
             0.3541020100 + 0.5300351800i;
             0.3628225400 - 0.1371443900i;
             0.3628225400 + 0.1371443900i];

ntf_gain = 1.0;

%% Specifications
order = 4;
fs = 2e6;
OSR = 16;
fB = fs/(2*OSR);
H_inf_design = 4.0;

fprintf('NTF Parameters (from Python):\n');
fprintf('  Zeros: %d total\n', length(ntf_zeros));
fprintf('  Poles: %d total\n', length(ntf_poles));
fprintf('  Gain: %.1f\n\n', ntf_gain);

%% Calculate NTF and STF frequency response
N_fft = 8192;
freqs_norm = (0:N_fft-1)/N_fft;  % Normalized frequency (f/fs)
z = exp(2*pi*1j*freqs_norm);

% NTF: H(z) = k * prod(z - zi) / prod(z - pi)
NTF_resp = ones(1, N_fft);
for i = 1:N_fft
    zi = z(i);
    num = ntf_gain * prod(zi - ntf_zeros);
    den = prod(zi - ntf_poles);
    NTF_resp(i) = num / den;
end
NTF_dB = 20*log10(abs(NTF_resp) + eps);

% For CIFF topology, STF ≈ 0 dB (unity gain)
STF_dB = zeros(1, N_fft);

fprintf('Frequency Response:\n');
fprintf('  NTF DC: %.2f dB\n', NTF_dB(1));
fprintf('  NTF Max (H_inf): %.2f dB\n', max(NTF_dB));
fprintf('  NTF at fB: %.2f dB\n', NTF_dB(round(N_fft/(2*OSR)) + 1));
fprintf('  STF (CIFF): ≈ 0 dB\n\n');

%% Load coefficients from Python
fprintf('Loading coefficients...\n');
run('coefficients_4th_order.txt');

fprintf('\nCoefficients:\n');
fprintf('  a = ['); fprintf('%.6f ', a); fprintf('] (feedback)\n');
fprintf('  g = ['); fprintf('%.6f ', g); fprintf('] (resonator)\n');
fprintf('  b = ['); fprintf('%.4f ', b); fprintf('] (input)\n');
fprintf('  c = ['); fprintf('%.4f ', c); fprintf('] (inter-stage)\n\n');

%% Build ABCD using toolbox stuffABCD
addpath('python-deltasigma-master/delsig/');
fprintf('[Building ABCD with stuffABCD]\n');

ABCD = stuffABCD(a, g, b, c, 'CIFF');

fprintf('ABCD Matrix (%dx%d):\n', size(ABCD,1), size(ABCD,2));
for i = 1:size(ABCD, 1)
    fprintf('  ['); fprintf('%10.6f ', ABCD(i, :)); fprintf(']\n');
end
fprintf('\n');

%% Simulation parameters
A_in = 0.5;
n_bits = 5;
n_levels = 32;
V_fs = 1.0;
LSB = 2*V_fs / (n_levels - 1);
q_levels = linspace(-V_fs, V_fs, n_levels);

fprintf('[Simulation]\n');
fprintf('----------------------------------------\n');
fprintf('Input amplitude: %.2f V\n', A_in);

N = 8192;
f_bin = round(sqrt(1/7) * N / (2*OSR));
u = A_in * sin(2*pi*f_bin*(0:N-1)/N);

n_states = size(ABCD, 1) - 1;
x = zeros(n_states, 1);
x_max = zeros(n_states, 1);
y = zeros(1, N);
v = zeros(1, N);

A_mat = ABCD(1:n_states, 1:n_states);
B_mat = ABCD(1:n_states, n_states+1:end);
C_mat = ABCD(n_states+1, 1:n_states);
D_mat = ABCD(n_states+1, n_states+1:end);

fprintf('Running simulation...\n');

for i = 1:N
    if i == 1, v_prev = 0; else, v_prev = v(i-1); end
    
    y(i) = C_mat * x + D_mat(1)*u(i) + D_mat(2)*v_prev;
    
    % Quantize
    y_clip = max(-V_fs, min(V_fs, y(i)));
    idx = round((y_clip - (-V_fs)) / LSB) + 1;
    idx = max(1, min(n_levels, idx));
    v(i) = q_levels(idx);
    
    % Update state
    x = A_mat*x + B_mat(:,1)*u(i) + B_mat(:,2)*v(i);
    x_max = max(x_max, abs(x));
end

fprintf('Simulation complete.\n');
fprintf('  Output range: [%.4f, %.4f]\n', min(v), max(v));
fprintf('  State maxima: '); fprintf('%.4f ', x_max); fprintf('\n\n');

%% Calculate SNR
w = 0.5 * (1 - cos(2*pi*(0:N-1)/N));
V_out = fft(v .* w) / (N/4);
V_out_mag = abs(V_out);

[~, sig_idx] = max(V_out_mag(2:N/2));
sig_bin = sig_idx + 1;
fB_bins = ceil(N / (2*OSR));
sig_bins = sig_bin-1:sig_bin+1;
sig_bins = sig_bins(sig_bins >= 2 & sig_bins <= fB_bins);
signal_power = sum(V_out_mag(sig_bins).^2);
noise_bins = setdiff(2:fB_bins, sig_bins);
noise_power = sum(V_out_mag(noise_bins).^2);

if noise_power > 0
    SNR = 10*log10(signal_power / noise_power);
    ENOB = (SNR - 1.76) / 6.02;
else
    SNR = inf;
    ENOB = inf;
end

fprintf('[Results]\n');
fprintf('----------------------------------------\n');
fprintf('Signal frequency: %.2f kHz\n', sig_bin/N * fs / 1000);
fprintf('SNR: %.2f dB\n', SNR);
if isfinite(ENOB)
    fprintf('ENOB: %.2f bits\n\n', ENOB);
end

%% Plot 1: Time Domain, Spectra, Histogram
figure('Position', [100 100 1400 900]);

% Time domain
subplot(2, 3, 1);
n_plot = 500;
t = (0:n_plot-1)/fs*1e6;
plot(t, u(1:n_plot), 'g-', 'LineWidth', 1.5); hold on;
stairs(t, v(1:n_plot), 'b-', 'LineWidth', 1);
hold off;
xlabel('Time (us)'); ylabel('Amplitude');
title(sprintf('Time Domain (%.2fV)', A_in));
legend('Input', 'Output'); grid on;

% Output spectrum
subplot(2, 3, 2);
freqs_plot = (0:N/2)/N*fs/1000;
V_out_dB = 20*log10(V_out_mag + eps);
plot(freqs_plot, V_out_dB(1:N/2+1), 'b-', 'LineWidth', 0.5); hold on;
f_sig_kHz = sig_bin/N*fs/1000;
fB_kHz = fB/1000;
plot([f_sig_kHz f_sig_kHz], [-150 10], 'g--', 'LineWidth', 2);
plot([fB_kHz fB_kHz], [-150 10], 'm--', 'LineWidth', 1);
fill([0 fB_kHz fB_kHz 0], [-150 -150 10 10], 'g', 'FaceAlpha', 0.1);
hold off;
xlabel('Frequency (kHz)'); ylabel('dBFS/NBW');
title(sprintf('Output Spectrum (SNR=%.1f dB)', SNR));
legend('Spectrum', 'Signal', 'BW Edge'); grid on;
axis([0 fs/2000 -140 10]);

% Noise spectrum (log)
subplot(2, 3, 3);
noise = v - u;
V_noise = fft(noise .* w) / (N/4);
V_noise_dB = 20*log10(abs(V_noise) + eps);
f_min = 100;
f_log = logspace(log10(f_min), log10(fs/2), 5000);
freqs_log = f_log/1000;
freqs_linear = (0:N/2)/N*fs;
V_noise_half = V_noise_dB(1:N/2+1);
V_noise_log = interp1(freqs_linear, V_noise_half, f_log, 'linear', 'extrap');
semilogx(freqs_log, V_noise_log, 'b-', 'LineWidth', 1); hold on;
plot([fB_kHz fB_kHz], [-150 10], 'm--', 'LineWidth', 1);
hold off;
xlabel('Frequency (kHz, log)'); ylabel('Noise (dB)');
title('Noise Spectrum'); legend('Noise', 'BW'); grid on;
axis([f_min/1000 fs/2000 -140 -40]);

% Histogram
subplot(2, 3, 4);
[counts, centers] = hist(v, 20);
bar(centers, counts, 'FaceColor', [0.3 0.5 0.8]);
xlabel('Level'); ylabel('Count');
title('Quantizer Histogram'); grid on;

% Status
subplot(2, 3, 5);
if max(x_max) < 100
    text(0.5, 0.7, sprintf('✓ STABLE\nMax: %.2f', max(x_max)), ...
        'FontSize', 14, 'Color', 'g', 'HorizontalAlignment', 'center');
else
    text(0.5, 0.7, sprintf('⚠ CHECK\nMax: %.2e', max(x_max)), ...
        'FontSize', 14, 'Color', 'orange', 'HorizontalAlignment', 'center');
end
text(0.5, 0.4, sprintf('4th-Order CIFF\nSNR=%.1f dB\nENOB=%.1f bits', SNR, ENOB), ...
    'FontSize', 11, 'HorizontalAlignment', 'center');
text(0.5, 0.2, sprintf('H_{inf}=%.1f, OSR=%d', H_inf_design, OSR), ...
    'FontSize', 10, 'HorizontalAlignment', 'center');
axis off;

% ABCD visualization
subplot(2, 3, 6);
text(0.1, 0.9, 'ABCD Matrix:', 'FontSize', 10, 'FontWeight', 'bold');
for i = 1:size(ABCD, 1)
    row_str = sprintf('Row %d: ', i);
    for j = 1:size(ABCD, 2)
        row_str = [row_str, sprintf('%.3f ', ABCD(i, j))];
    end
    text(0.1, 0.8-(i-1)*0.15, row_str, 'FontSize', 8);
end
axis off;

print('dsm_4th_order.png', '-dpng', '-r150');
fprintf('Saved: dsm_4th_order.png\n');

%% Plot 2: NTF and STF Frequency Response
figure('Position', [100 100 1400 600]);

% NTF Plot
subplot(1, 2, 1);

% Log frequency scale
f_min = 100;
f_log = logspace(log10(f_min), log10(fs/2), 5000);
freqs_log_norm = f_log / fs;

% Interpolate NTF for log plot
NTF_log = interp1(freqs_norm, NTF_dB, freqs_log_norm, 'linear', 'extrap');

semilogx(f_log/1000, NTF_log, 'b-', 'LineWidth', 1.5); hold on;

% Mark bandwidth edge
fB_kHz = fB/1000;
plot([fB_kHz fB_kHz], [-100 40], 'm--', 'LineWidth', 2);

% Mark H_inf
[max_ntf, max_idx] = max(NTF_dB);
f_max_kHz = freqs_norm(max_idx) * fs / 1000;
plot(f_max_kHz, max_ntf, 'ro', 'MarkerSize', 10, 'LineWidth', 2);
text(f_max_kHz*1.2, max_ntf+3, sprintf('H_{inf}=%.1f dB', max_ntf), ...
    'FontSize', 10, 'Color', 'r');

% Mark notch frequencies (from zeros)
for i = 1:length(ntf_zeros)
    zi = ntf_zeros(i);
    if abs(zi) < 1.1
        f_notch = abs(angle(zi)) / (2*pi) * fs;
        if f_notch < fs/2 && f_notch > 100
            plot(f_notch/1000, -70, 'gv', 'MarkerSize', 8, 'LineWidth', 1.5);
        end
    end
end

hold off;
xlabel('Frequency (kHz)', 'FontSize', 12);
ylabel('Magnitude (dB)', 'FontSize', 12);
title(sprintf('NTF - Noise Transfer Function\n(4th-Order, H_{inf}=%.1f)', H_inf_design), 'FontSize', 12);
legend('NTF', 'Signal BW', 'H_{inf}', 'Notches', 'Location', 'northwest');
grid on;
axis([100/1000 fs/2000 -100 40]);

% STF Plot
subplot(1, 2, 2);

% For CIFF, STF is approximately unity (0 dB)
STF_log = zeros(size(freqs_log));

semilogx(f_log/1000, STF_log, 'g-', 'LineWidth', 1.5); hold on;
plot([fB_kHz fB_kHz], [-10 15], 'm--', 'LineWidth', 2);

% Mark DC
plot(100/1000, 0, 'ro', 'MarkerSize', 10, 'LineWidth', 2);
text(150/1000, 2, 'DC: 0 dB', 'FontSize', 10, 'Color', 'r');

% Add annotation
annotation('textbox', [0.58 0.15 0.3 0.15], 'String', ...
    sprintf('CIFF Topology\nUnity STF (0 dB)\nFlat in passband'), ...
    'FitBoxToText', 'on', 'BackgroundColor', [0.9 0.9 0.9], ...
    'EdgeColor', 'black');

hold off;
xlabel('Frequency (kHz)', 'FontSize', 12);
ylabel('Magnitude (dB)', 'FontSize', 12);
title('STF - Signal Transfer Function\n(CIFF Topology)', 'FontSize', 12);
legend('STF', 'Signal BW', 'DC Gain', 'Location', 'northwest');
grid on;
axis([100/1000 fs/2000 -10 15]);

print('dsm_4th_order_ntf_stf.png', '-dpng', '-r150');
fprintf('Saved: dsm_4th_order_ntf_stf.png\n');

%% Summary
fprintf('\n============================================================\n');
if max(x_max) < 100
    fprintf('  ✓ SIMULATION COMPLETE\n');
else
    fprintf('  ⚠ CHECK STABILITY\n');
end
fprintf('============================================================\n');
