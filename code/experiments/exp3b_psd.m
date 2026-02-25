% exp3b_psd.m - 实验3b：波形功率谱密度对比（复现论文 Fig.8/9）
%
% 对比 PRD 波形、文献[17]波形和目标 PSD 的频域特性：
%   高SINR：PRD 波形 PSD 覆盖目标 PSD 的整体形状（宽频）
%   低SINR：两种方法接近一致（窄带，集中于最强频点）
%
% 目标1：sigma2 ∈ {-10, -5, 0, 5} dB  →  Fig.8
% 目标2：sigma2 ∈ {-10, 0, 10, 20} dB →  Fig.9  （注：目标2低SINR sigma2较大）

clear; clc; close all;
script_dir = fileparts(mfilename('fullpath'));
code_dir   = fileparts(script_dir);
addpath(fullfile(code_dir, 'utils'));
addpath(fullfile(code_dir, 'algorithms'));

gamma0_dB = 15;
gamma0    = 10^(gamma0_dB / 10);
% rho 在循环中设为 1/L (CMC)，与论文 Section V.C 一致
Nfft      = 512;

target_configs = {
    1, [-10, -5, 0, 5];
    2, [-10, 0, 10, 20]
};
fig_letters = {{'(a)','(b)','(c)','(d)'}, {'(a)','(b)','(c)','(d)'}};
fig_nums    = [8, 9];

for tc = 1:2
    target_id     = target_configs{tc, 1};
    sigma2_dB_vec = target_configs{tc, 2};

    p = get_target_params(target_id);
    L = p.L;  Q = p.Q;  t0 = p.t0;  Rt = p.Rt;
    rho = 1/L;  % CMC 恒模约束

    % 计算目标 PSD（基于 t0 的频谱作为参考）
    t0_padded = [t0; zeros(Nfft - Q, 1)];
    target_psd_dB = 20 * log10(abs(fftshift(fft(t0_padded, Nfft))) + 1e-12);
    target_psd_dB = target_psd_dB - max(target_psd_dB);  % 归一化到 0 dB

    freq_axis = linspace(0, 1, Nfft);   % 归一化频率

    figure('Position', [100, 100, 1200, 900]);

    for si = 1:length(sigma2_dB_vec)
        sigma2 = 10^(sigma2_dB_vec(si) / 10);
        Rn_inv = (1/sigma2) * eye(L + Q - 1);

        fprintf('[目标%d] sigma2=%ddB，设计波形...\n', target_id, sigma2_dB_vec(si));

        % PRD 方法
        [s_prd, ~, ~, ~, ~, ~] = PRD_design(t0, Rt, Rn_inv, rho, L, Q, gamma0);

        % 文献[17]方法
        [s_ref, ~, ~, ~] = baseline_metric17(t0, Rt, Rn_inv, rho, L, Q, gamma0);

        % 计算 PSD（FFT，dB，归一化）
        s_prd_pad = [s_prd; zeros(Nfft - L, 1)];
        s_ref_pad = [s_ref; zeros(Nfft - L, 1)];

        psd_prd = 20 * log10(abs(fftshift(fft(s_prd_pad, Nfft))) + 1e-12);
        psd_ref = 20 * log10(abs(fftshift(fft(s_ref_pad, Nfft))) + 1e-12);

        % 各自归一化
        psd_prd = psd_prd - max(psd_prd);
        psd_ref = psd_ref - max(psd_ref);

        subplot(2, 2, si);
        hold on; grid on;
        plot(freq_axis, target_psd_dB, 'k-',  'LineWidth', 2,   'DisplayName', 'Target PSD');
        plot(freq_axis, psd_prd,       'b-',  'LineWidth', 1.5, 'DisplayName', 'Our metric (PRD)');
        plot(freq_axis, psd_ref,       'r-',  'LineWidth', 1.5, 'DisplayName', 'Metric in [17]');

        ylim([-40, 5]);
        legend('show', 'FontSize', 8, 'Location', 'southwest');
        xlabel('Normalized Frequency');
        ylabel('Power/dB');
        title(sprintf('%s \\sigma^2 = %d dB', fig_letters{tc}{si}, sigma2_dB_vec(si)));
    end

    sgtitle(sprintf('Fig.%d: Waveform PSD for Target %d', fig_nums(tc), target_id), 'FontSize', 13);
    fname = sprintf('result_fig%d_psd_target%d.png', fig_nums(tc), target_id);
    saveas(gcf, fname);
    fprintf('目标%d PSD图已保存为 %s\n', target_id, fname);
end
