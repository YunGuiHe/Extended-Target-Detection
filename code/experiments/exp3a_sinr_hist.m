% exp3a_sinr_hist.m - 实验3a：SINR 分布直方图对比（复现论文 Fig.6/7）
%
% 对比 PRD 方法和文献[17]方法在随机 TIR 下的 SINR 分布：
%   1000 次蒙特卡罗采样 t~CN(t0, Rt)，计算 gamma(s) = t^H*M*t
%   绘制双色柱状直方图，竖虚线标出目标 SINR 阈值 gamma0
%
% 目标1：sigma2 ∈ {-10, -5, 0, 5} dB
% 目标2：sigma2 ∈ {-10, 0, 10, 20} dB

clear; clc; close all;
script_dir = fileparts(mfilename('fullpath'));
code_dir   = fileparts(script_dir);
addpath(fullfile(code_dir, 'utils'));
addpath(fullfile(code_dir, 'algorithms'));

rng(42);
gamma0_dB = 15;
gamma0    = 10^(gamma0_dB / 10);
N_MC      = 1000;   % 蒙特卡罗次数
% 论文 Section V.C: "we set rho = 1/L (CMC)"
% 恒模约束（Constant Modulus Constraint）是Fig.6/7的正确参数

target_configs = {
    1, [-10, -5, 0, 5],   {'(a)', '(b)', '(c)', '(d)'};
    2, [-10, 0, 10, 20],  {'(a)', '(b)', '(c)', '(d)'}
};

for tc = 1:2
    target_id     = target_configs{tc, 1};
    sigma2_dB_vec = target_configs{tc, 2};
    sub_labels    = target_configs{tc, 3};

    p = get_target_params(target_id);
    L = p.L;  Q = p.Q;  t0 = p.t0;  Rt = p.Rt;
    rho = 1/L;  % CMC 恒模约束（论文 Section V.C）

    % Cholesky 分解 Rt 用于生成随机 TIR 样本
    Rt_chol = chol(Rt + 1e-9*eye(Q), 'lower');

    figure('Position', [100, 100, 1200, 900]);

    for si = 1:length(sigma2_dB_vec)
        sigma2 = 10^(sigma2_dB_vec(si) / 10);
        Rn_inv = (1/sigma2) * eye(L + Q - 1);

        fprintf('[目标%d] sigma2=%ddB，设计波形...\n', target_id, sigma2_dB_vec(si));

        % PRD 方法设计波形
        [s_prd, Rs_prd, f_Rs_prd, f_s_prd, ~, sinr_case] = PRD_design(t0, Rt, Rn_inv, rho, L, Q, gamma0);

        % 文献[17]方法设计波形
        [s_ref, Rs_ref, f_Rs_ref, f_s_ref] = baseline_metric17(t0, Rt, Rn_inv, rho, L, Q, gamma0);

        % 诊断：检查两种方法的波形差异
        s_diff = norm(abs(s_prd) - abs(s_ref)) / max(norm(s_prd), 1e-12);
        Rs_diff = norm(Rs_prd - Rs_ref, 'fro') / max(norm(Rs_prd, 'fro'), 1e-12);
        fprintf('  诊断：||Rs_prd-Rs_ref||/||Rs_prd|| = %.6f\n', Rs_diff);
        fprintf('  诊断：||s_prd|-|s_ref||/||s_prd|| = %.6f\n', s_diff);
        fprintf('  诊断：f_prd(Rs)=%.4f, f_ref(Rs)=%.4f\n', f_Rs_prd, f_Rs_ref);
        fprintf('  诊断：f_prd(s)=%.4f, f_ref(s)=%.4f\n', f_s_prd, f_s_ref);
        fprintf('  诊断：SINR case = %s\n', sinr_case);

        % 构造波形矩阵 S 并计算 M
        coeffs_prd = precompute_coeffs(t0, Rt, Rn_inv, Q, L);
        F = coeffs_prd.F;
        M_prd = compute_M(s_prd * s_prd', F, Q);
        M_ref = compute_M(s_ref * s_ref', F, Q);

        % 诊断：检查 M 矩阵差异
        M_diff = norm(M_prd - M_ref, 'fro') / max(norm(M_prd, 'fro'), 1e-12);
        fprintf('  诊断：||M_prd-M_ref||/||M_prd|| = %.6f\n', M_diff);

        % 计算两种方法的 E[gamma] 和 D[gamma]
        [E_prd, V_prd, ~] = compute_EV_from_M(M_prd, t0, Rt, gamma0);
        [E_ref, V_ref, ~] = compute_EV_from_M(M_ref, t0, Rt, gamma0);
        fprintf('  PRD:  E[gamma]=%.2f (%.1f dB), D[gamma]=%.2f (%.1f dB)\n', ...
            E_prd, 10*log10(max(E_prd,1e-10)), sqrt(max(V_prd,0)), 10*log10(max(sqrt(max(V_prd,0)),1e-10)));
        fprintf('  [17]: E[gamma]=%.2f (%.1f dB), D[gamma]=%.2f (%.1f dB)\n', ...
            E_ref, 10*log10(max(E_ref,1e-10)), sqrt(max(V_ref,0)), 10*log10(max(sqrt(max(V_ref,0)),1e-10)));

        % 蒙特卡罗：生成 1000 个随机 TIR 样本
        gamma_prd = zeros(N_MC, 1);
        gamma_ref = zeros(N_MC, 1);
        for mc = 1:N_MC
            t_rand = t0 + Rt_chol * (randn(Q, 1) + 1i * randn(Q, 1)) / sqrt(2);
            gamma_prd(mc) = real(t_rand' * M_prd * t_rand);
            gamma_ref(mc) = real(t_rand' * M_ref * t_rand);
        end

        % 转为 dB
        gamma_prd_dB = 10 * log10(max(gamma_prd, 1e-10));
        gamma_ref_dB = 10 * log10(max(gamma_ref, 1e-10));

        subplot(2, 2, si);
        hold on;
        all_vals = [gamma_prd_dB; gamma_ref_dB];
        bin_edges = linspace(min(all_vals) - 1, max(all_vals) + 1, 31);

        % 使用 histcounts 分别计算两个分布的频数
        [counts_prd, ~] = histcounts(gamma_prd_dB, bin_edges);
        [counts_ref, ~] = histcounts(gamma_ref_dB, bin_edges);

        % 将双色柱状图并排显示（Grouped Bar Chart）以还原论文 Fig.6/7 的展现形式
        % 为了让柱子对齐到 bin 的中心，计算 bin_centers
        bin_centers = bin_edges(1:end-1) + diff(bin_edges)/2;
        
        % 组合频数矩阵，并用 bar 函数绘制 'grouped' 模式
        counts_matrix = [counts_prd(:), counts_ref(:)]; 
        b = bar(bin_centers, counts_matrix, 'grouped', 'BarWidth', 1);
        
        % 设置颜色以匹配论文 (蓝色代表 Our metric (PRD)，红色代表 Metric in [17])
        b(1).FaceColor = 'b';
        b(2).FaceColor = 'r';

        % 标出 gamma0
        xline(gamma0_dB, 'k--', '\gamma_0', 'LineWidth', 2, 'LabelOrientation', 'horizontal');

        % 手动添加图例
        legend([b(1), b(2)], {'Our metric (PRD)', 'Metric in [17]'}, 'FontSize', 9, 'Location', 'northwest');
        xlabel('SINR/dB');
        ylabel('Frequency');
        title(sprintf('%s \\sigma^2 = %d dB', sub_labels{si}, sigma2_dB_vec(si)));
        grid on;
    end

    sgtitle(sprintf('Fig.%d: SINR Histogram for Target %d', 5 + tc, target_id), 'FontSize', 13);
    fname = sprintf('result_fig%d_sinr_hist_target%d.png', 5 + tc, target_id);
    saveas(gcf, fname);
    fprintf('目标%d SINR直方图已保存为 %s\n', target_id, fname);
end
