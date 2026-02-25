% exp4_Psd_curve.m - 实验5：稳定检测概率曲线（复现论文 Fig.10）
%
% 对比 PRD 和 [17] 两种方法的稳定检测概率 Psd：
%   对每个 sigma2，设计两种波形，蒙特卡洛生成 TIR 样本，
%   统计 gamma >= gamma0 的比例作为 Psd
%
% 目标1：sigma2 从 -10 dB 到 8 dB（细粒度扫描）
% 目标2：sigma2 从 -10 dB 到 25 dB

clear; clc; close all;
script_dir = fileparts(mfilename('fullpath'));
code_dir   = fileparts(script_dir);
addpath(fullfile(code_dir, 'utils'));
addpath(fullfile(code_dir, 'algorithms'));

rng(42);
gamma0_dB = 15;
gamma0    = 10^(gamma0_dB / 10);
N_MC      = 1000;   % 蒙特卡洛次数（论文 Section V.C 使用 1000）

% 目标配置：扫描的 sigma2 范围
target_configs = {
    1, -10:2:8;     % 目标1：sigma2 从 -10 到 8 dB，步长 2 dB
    2, -10:3:25     % 目标2：sigma2 从 -10 到 25 dB，步长 3 dB
};

figure('Position', [100, 100, 1000, 450]);

for tc = 1:2
    target_id     = target_configs{tc, 1};
    sigma2_dB_vec = target_configs{tc, 2};
    n_sigma       = length(sigma2_dB_vec);

    p = get_target_params(target_id);
    L = p.L;  Q = p.Q;  t0 = p.t0;  Rt = p.Rt;
    rho = 1/L;  % CMC 恒模约束

    % Cholesky 分解 Rt 用于生成随机 TIR 样本
    Rt_chol = chol(Rt + 1e-9*eye(Q), 'lower');

    Psd_prd = zeros(1, n_sigma);
    Psd_ref = zeros(1, n_sigma);

    for si = 1:n_sigma
        sigma2 = 10^(sigma2_dB_vec(si) / 10);
        Rn_inv = (1/sigma2) * eye(L + Q - 1);

        fprintf('[目标%d] sigma2=%ddB (%d/%d)，设计波形...\n', ...
                target_id, sigma2_dB_vec(si), si, n_sigma);

        % PRD 方法设计波形
        [s_prd, ~, ~, ~, ~, sinr_case] = PRD_design(t0, Rt, Rn_inv, rho, L, Q, gamma0);

        % 文献[17]方法设计波形
        [s_ref, ~, ~, ~] = baseline_metric17(t0, Rt, Rn_inv, rho, L, Q, gamma0);

        % 计算 M 矩阵
        coeffs_tmp = precompute_coeffs(t0, Rt, Rn_inv, Q, L);
        F = coeffs_tmp.F;
        M_prd = compute_M(s_prd * s_prd', F, Q);
        M_ref = compute_M(s_ref * s_ref', F, Q);

        % 蒙特卡洛仿真：统计 gamma >= gamma0 的次数
        count_prd = 0;
        count_ref = 0;
        for mc = 1:N_MC
            t_rand = t0 + Rt_chol * (randn(Q, 1) + 1i * randn(Q, 1)) / sqrt(2);
            gamma_prd = real(t_rand' * M_prd * t_rand);
            gamma_ref = real(t_rand' * M_ref * t_rand);
            if gamma_prd >= gamma0
                count_prd = count_prd + 1;
            end
            if gamma_ref >= gamma0
                count_ref = count_ref + 1;
            end
        end

        Psd_prd(si) = count_prd / N_MC;
        Psd_ref(si) = count_ref / N_MC;

        fprintf('  SINR case = %s, Psd_prd = %.3f, Psd_ref = %.3f\n', ...
                sinr_case, Psd_prd(si), Psd_ref(si));
    end

    % 绘图
    subplot(1, 2, tc);
    hold on; grid on;
    plot(sigma2_dB_vec, Psd_prd, '-bs', 'LineWidth', 1.5, 'MarkerSize', 6, ...
         'MarkerFaceColor', 'b', 'DisplayName', 'Our metric (PRD)');
    plot(sigma2_dB_vec, Psd_ref, '-ro', 'LineWidth', 1.5, 'MarkerSize', 6, ...
         'MarkerFaceColor', 'r', 'DisplayName', 'Metric in [17]');

    xlabel('\sigma^2 / dB', 'FontSize', 12);
    ylabel('P_{sd}', 'FontSize', 12);
    title(sprintf('(%s) Target %d', char('a' + tc - 1), target_id), 'FontSize', 13);
    legend('show', 'FontSize', 10, 'Location', 'northwest');
    ylim([-0.05, 1.05]);
    set(gca, 'FontSize', 11);
end

sgtitle('Fig.10: P_{sd} versus \sigma^2', 'FontSize', 14);
saveas(gcf, 'result_fig10_Psd.png');
fprintf('\n实验5完成，结果已保存为 result_fig10_Psd.png\n');
