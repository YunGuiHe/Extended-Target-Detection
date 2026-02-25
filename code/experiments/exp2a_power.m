% exp2a_power.m - 实验2a：波形功率曲线（复现论文 Fig.4）
%
% 验证合成波形每个采样点的功率 |s*(i)|^2 均不超过 PAR 上限 rho
%
% 参数（对应论文 Fig.4）：
%   目标1：rho = 1/4 (sigma2=-3dB) 和 rho = 1/20 (sigma2=-3dB)
%   目标2：rho = 1/30 (sigma2=-6dB) 和 rho = 1/L(CMC) (sigma2=-6dB)

clear; clc; close all;
script_dir = fileparts(mfilename('fullpath'));
code_dir   = fileparts(script_dir);
addpath(fullfile(code_dir, 'utils'));
addpath(fullfile(code_dir, 'algorithms'));

gamma0_dB = 15;
gamma0    = 10^(gamma0_dB / 10);

figure('Position', [100, 100, 900, 700]);

%% 子图(a)：目标1
subplot(2, 1, 1);
p = get_target_params(1);
L = p.L;  Q = p.Q;  t0 = p.t0;  Rt = p.Rt;
sigma2 = 10^(-3/10);
Rn_inv = (1/sigma2) * eye(L + Q - 1);

configs = {1/4, '-bs', 'PMR = 1/4,  \sigma^2 = -3dB'; ...
           1/20, '-ro', '\rho = 1/20, \sigma^2 = -3dB'};

hold on; grid on;
for ci = 1:size(configs, 1)
    rho    = configs{ci, 1};
    lstyle = configs{ci, 2};
    lname  = configs{ci, 3};

    fprintf('[目标1] rho=%.4f, sigma2=-3dB...\n', rho);
    [s_star, ~, ~, ~, ~, ~] = PRD_design(t0, Rt, Rn_inv, rho, L, Q, gamma0);
    power_wf = abs(s_star).^2;
    plot(1:L, power_wf, lstyle, 'LineWidth', 1.5, 'MarkerSize', 5, ...
         'DisplayName', lname);

    % 验证 PAR 约束
    fprintf('  最大功率 = %.5f，PAR上限 = %.5f，满足约束：%d\n', ...
            max(power_wf), rho, max(power_wf) <= rho + 1e-6);
end
legend('show', 'FontSize', 10);
xlabel('Pulse samples');
ylabel('Power');
title('(a) Power of Designed Waveform: Target 1');
xlim([1, L]);

%% 子图(b)：目标2
subplot(2, 1, 2);
p = get_target_params(2);
L = p.L;  Q = p.Q;  t0 = p.t0;  Rt = p.Rt;
sigma2 = 10^(-6/10);
Rn_inv = (1/sigma2) * eye(L + Q - 1);

configs2 = {1/30,  '-bs', '\rho = 1/30, \sigma^2 = -6dB'; ...
            1/L,   '-ro', '\rho = 1/L (CMC), \sigma^2 = -6dB'};

hold on; grid on;
for ci = 1:size(configs2, 1)
    rho    = configs2{ci, 1};
    lstyle = configs2{ci, 2};
    lname  = configs2{ci, 3};

    fprintf('[目标2] rho=%.4f, sigma2=-6dB...\n', rho);
    [s_star, ~, ~, ~, ~, ~] = PRD_design(t0, Rt, Rn_inv, rho, L, Q, gamma0);
    power_wf = abs(s_star).^2;
    plot(1:L, power_wf, lstyle, 'LineWidth', 1.5, 'MarkerSize', 5, ...
         'DisplayName', lname);

    fprintf('  最大功率 = %.5f，PAR上限 = %.5f，满足约束：%d\n', ...
            max(power_wf), rho, max(power_wf) <= rho + 1e-6);
end
legend('show', 'FontSize', 10);
xlabel('Pulse samples');
ylabel('Power');
title('(b) Power of Designed Waveform: Target 2');
xlim([1, L]);

sgtitle('Fig.4: Power Curves of Designed Waveforms', 'FontSize', 13);
saveas(gcf, 'result_fig4_power.png');
fprintf('实验2a完成，结果已保存为 result_fig4_power.png\n');
