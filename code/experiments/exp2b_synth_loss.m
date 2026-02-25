% exp2b_synth_loss.m - 实验2b：合成损失（复现论文 Fig.5）
%
% 扫描 rho，比较 f(Rs)（虚线/矩阵最优值）与 f(s*)（实线/合成波形值）
% 验证随机化合成损失 < 0.1
%
% 参数：
%   目标1：sigma2 = {-3, 3} dB
%   目标2：sigma2 = {-6, 10} dB

clear; clc; close all;
script_dir = fileparts(mfilename('fullpath'));
code_dir   = fileparts(script_dir);
addpath(fullfile(code_dir, 'utils'));
addpath(fullfile(code_dir, 'algorithms'));

gamma0_dB = 15;
gamma0    = 10^(gamma0_dB / 10);

figure('Position', [100, 100, 900, 700]);

target_configs = {
    1, [-3, 3],   linspace(1/60, 1, 15);     % 目标1，rho 扫描
    2, [-6, 10],  linspace(1/60, 1/1.5, 15)  % 目标2，rho 扫描
};

colors = {'b', 'r'};

for tc = 1:2
    target_id     = target_configs{tc, 1};
    sigma2_dB_vec = target_configs{tc, 2};
    rho_vec       = target_configs{tc, 3};

    p = get_target_params(target_id);
    L = p.L;  Q = p.Q;  t0 = p.t0;  Rt = p.Rt;

    subplot(2, 1, tc);
    hold on; grid on;

    for si = 1:length(sigma2_dB_vec)
        sigma2 = 10^(sigma2_dB_vec(si) / 10);
        Rn_inv = (1/sigma2) * eye(L + Q - 1);
        col    = colors{si};

        f_Rs_all = zeros(size(rho_vec));
        f_s_all  = zeros(size(rho_vec));

        for ri = 1:length(rho_vec)
            rho = rho_vec(ri);
            fprintf('[目标%d] sigma2=%ddB, rho=%.4f...\n', target_id, sigma2_dB_vec(si), rho);

            coeffs = precompute_coeffs(t0, Rt, Rn_inv, Q, L);
            [Rs0, f0, obj0] = solve_init_SDP(t0, Rt, rho, L, Q, gamma0, coeffs);

            if obj0 >= 0
                [Rs_opt, ~] = solve_high_SINR(t0, Rt, rho, L, Q, gamma0, coeffs, Rs0, f0);
            else
                [Rs_opt, ~] = solve_low_SINR(t0, Rt, rho, L, Q, gamma0, coeffs, Rs0, f0);
            end

            M_opt = compute_M(Rs_opt, coeffs.F, Q);
            [~, ~, f_Rs] = compute_EV_from_M(M_opt, t0, Rt, gamma0);
            f_Rs_all(ri) = f_Rs;

            [~, f_s] = synthesize_waveform(Rs_opt, t0, Rt, coeffs, Q, L, gamma0, rho);
            f_s_all(ri) = f_s;
        end

        % 虚线：f(Rs)，实线：f(s*)
        plot(rho_vec, f_Rs_all, ['--', col], 'LineWidth', 1.5, ...
             'DisplayName', sprintf('f(R_s), \\sigma^2=%ddB', sigma2_dB_vec(si)));
        plot(rho_vec, f_s_all,  ['-',  col], 'LineWidth', 1.5, 'Marker', 'o', 'MarkerSize', 5, ...
             'DisplayName', sprintf('f(s^*), \\sigma^2=%ddB', sigma2_dB_vec(si)));

        fprintf('  最大合成损失 = %.4f\n', max(abs(f_Rs_all - f_s_all)));
    end

    legend('show', 'FontSize', 9, 'Location', 'best');
    xlabel('\rho');
    ylabel('Objective value');
    title(sprintf('(%s) Synthesizing Losses: Target %d', char('a' + tc - 1), target_id));
end

sgtitle('Fig.5: Synthesizing Losses vs. PAR Constraint', 'FontSize', 13);
saveas(gcf, 'result_fig5_synth_loss.png');
fprintf('实验2b完成，结果已保存为 result_fig5_synth_loss.png\n');
