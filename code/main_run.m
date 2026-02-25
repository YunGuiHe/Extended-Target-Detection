% main_run.m - 主运行脚本
% 按顺序运行论文复现的前三组实验
%
% 依赖：CVX 工具箱（运行前请确保已安装并执行 cvx_setup）
%
% 目录结构：
%   code/
%   ├── utils/          (compute_J, get_target_params, precompute_coeffs, ...)
%   ├── algorithms/     (solve_init_SDP, solve_high_SINR, solve_low_SINR, ...)
%   └── experiments/    (exp1_convergence, exp2a_power, ...)

clear; clc; close all;

% 添加路径
addpath('utils', 'algorithms', 'experiments');

% 检查 CVX 是否安装
if ~exist('cvx_begin', 'file')
    error('未检测到 CVX 工具箱，请先运行 cvx_setup');
end

fprintf('==========================================\n');
fprintf('  PRD 波形设计论文复现实验\n');
fprintf('==========================================\n\n');

% 切换到 experiments 目录存放结果图片
orig_dir = pwd;
cd('experiments');

%% 实验1：收敛性验证（Fig.3）
% 验证算法在 5 次以内收敛
run_exp('exp1_convergence', '实验1：收敛性验证（Fig.3）');

%% 实验2a：波形功率曲线（Fig.4）
% 验证合成波形满足 PAR 约束
run_exp('exp2a_power', '实验2a：波形功率曲线（Fig.4）');

%% 实验2b：合成损失（Fig.5）
% 验证 SDR 松弛的合成损失 < 0.1
run_exp('exp2b_synth_loss', '实验2b：合成损失（Fig.5）');

%% 实验3a：SINR 分布直方图（Fig.6 / Fig.7）
% 对比 PRD 和 [17] 方法的 SINR 分布
run_exp('exp3a_sinr_hist', '实验3a：SINR直方图（Fig.6/7）');

%% 实验3b：波形功率谱密度（Fig.8 / Fig.9）
% 频域对比两种波形设计方法
run_exp('exp3b_psd', '实验3b：波形PSD（Fig.8/9）');

%% 实验4：稳定检测概率（Fig.10）
% 定量评估 PRD vs [17] 的检测概率优势
run_exp('exp4_Psd_curve', '实验4：稳定检测概率（Fig.10）');

cd(orig_dir);

fprintf('\n==========================================\n');
fprintf('  所有实验完成！结果图片保存在 experiments/ 目录\n');
fprintf('==========================================\n');

%% 辅助函数
function run_exp(exp_name, desc)
    fprintf('\n------------------------------------------\n');
    fprintf('运行：%s\n', desc);
    fprintf('------------------------------------------\n');
    t_start = tic;
    try
        run(exp_name);
        fprintf('[完成] 耗时 %.1f 秒\n', toc(t_start));
    catch ME
        fprintf('[错误] %s\n', ME.message);
        fprintf('  文件：%s，行：%d\n', ME.stack(1).file, ME.stack(1).line);
    end
end
