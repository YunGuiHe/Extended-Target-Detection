% exp1_convergence.m - 实验1：收敛性验证（复现论文 Fig.3）
%
% 对目标1和目标2，分别在高/低SINR情况下运行算法，
% 记录每次 Dinkelbach 辅助目标值，绘制收敛曲线（4 个子图）
%
% 对应论文参数：
%   高SINR：sigma2_dB = {-3, -6} dB，rho = {1/4, 1/L}
%   低SINR（目标1）：sigma2_dB = {3, 6} dB
%   低SINR（目标2）：sigma2_dB = {10, 13} dB

clear; clc; close all;
% 基于脚本所在目录动态计算路径（兼容任意当前工作目录）
script_dir = fileparts(mfilename('fullpath'));
code_dir   = fileparts(script_dir);           % code/ 目录
addpath(fullfile(code_dir, 'utils'));
addpath(fullfile(code_dir, 'algorithms'));

gamma0_dB = 15;
gamma0    = 10^(gamma0_dB / 10);

% ======== 高SINR 参数 ========
sigma2_high_dB = [-3, -6];
rho_vals_high  = {'1/4', '1/L'};

% ======== 低SINR 参数 ========
sigma2_low_T1_dB = [3, 6];    % 目标1
sigma2_low_T2_dB = [10, 13];  % 目标2

labels_high  = {'$\rho=1/4,\,\sigma^2=-3$dB', '$\rho=1/L\,(CMC),\,\sigma^2=-3$dB', ...
                '$\rho=1/4,\,\sigma^2=-6$dB', '$\rho=1/L\,(CMC),\,\sigma^2=-6$dB'};
line_styles  = {'-bs', '-ro', '--b^', '--r*'};

figure('Position', [100, 100, 1200, 900]);

for target_id = 1:2
    p = get_target_params(target_id);
    L = p.L;  Q = p.Q;  t0 = p.t0;  Rt = p.Rt;

    if target_id == 1
        sigma2_low_dB = sigma2_low_T1_dB;
        subplot_high  = (target_id - 1) * 2 + 1;   % 子图(a)
        subplot_low   = (target_id - 1) * 2 + 2;   % 子图(b)
    else
        sigma2_low_dB = sigma2_low_T2_dB;
        subplot_high  = (target_id - 1) * 2 + 1;   % 子图(c)
        subplot_low   = (target_id - 1) * 2 + 2;   % 子图(d)
    end

    % ---- 高SINR ----
    subplot(2, 2, subplot_high);
    hold on; grid on;
    legend_str = {};
    li = 1;
    for sig_idx = 1:length(sigma2_high_dB)
        sigma2 = 10^(sigma2_high_dB(sig_idx) / 10);
        Rn_inv = (1/sigma2) * eye(L + Q - 1);

        for ri = 1:2
            if ri == 1
                rho = 1/4;
            else
                rho = 1/L;
            end

            fprintf('[目标%d 高SINR] sigma2=%ddB, rho=%s...\n', ...
                    target_id, sigma2_high_dB(sig_idx), rho_vals_high{ri});
            [~, ~, ~, ~, f_vals, ~] = PRD_design(t0, Rt, Rn_inv, rho, L, Q, gamma0);

            iters = 1:length(f_vals);
            plot(iters, f_vals, line_styles{li}, 'LineWidth', 1.5);
            legend_str{end+1} = labels_high{(sig_idx-1)*2 + ri};  %#ok<SAGROW>
            li = li + 1;
        end
    end
    legend(legend_str, 'Interpreter', 'latex', 'FontSize', 8, 'Location', 'northeast');
    xlabel('Number of Iterations');
    ylabel('Objective Value');
    title(sprintf('((%s)) Target %d, High SINR', char('a' + (target_id-1)*2), target_id));

    % ---- 低SINR ----
    subplot(2, 2, subplot_low);
    hold on; grid on;
    legend_str = {};
    li = 1;
    for sig_idx = 1:length(sigma2_low_dB)
        sigma2 = 10^(sigma2_low_dB(sig_idx) / 10);
        Rn_inv = (1/sigma2) * eye(L + Q - 1);

        for ri = 1:2
            if ri == 1
                rho = 1/4;
            else
                rho = 1/L;
            end

            fprintf('[目标%d 低SINR] sigma2=%ddB, rho=%s...\n', ...
                    target_id, sigma2_low_dB(sig_idx), rho_vals_high{ri});
            [~, ~, ~, ~, f_vals, ~] = PRD_design(t0, Rt, Rn_inv, rho, L, Q, gamma0);

            label = sprintf('$\\rho=%s,\\,\\sigma^2=%d$dB', rho_vals_high{ri}, sigma2_low_dB(sig_idx));
            plot(1:length(f_vals), f_vals, line_styles{li}, 'LineWidth', 1.5);
            legend_str{end+1} = label;  %#ok<SAGROW>
            li = li + 1;
        end
    end
    legend(legend_str, 'Interpreter', 'latex', 'FontSize', 8, 'Location', 'northeast');
    xlabel('Number of Iterations');
    ylabel('Objective Value');
    title(sprintf('((%s)) Target %d, Low SINR', char('a' + (target_id-1)*2 + 1), target_id));
    if target_id == 1
        xlim([1, 4]);   % 论文 Fig.3(b) 横轴范围
    else
        xlim([1, 3]);   % 论文 Fig.3(d) 横轴范围
    end
end

sgtitle('Fig.3: Iteration Curves - Convergence Verification', 'FontSize', 14);
saveas(gcf, 'result_fig3_convergence.png');
fprintf('实验1完成，结果已保存为 result_fig3_convergence.png\n');
