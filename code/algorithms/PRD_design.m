function [s_star, Rs_opt, f_Rs, f_s, f_vals, sinr_case] = PRD_design(t0, Rt, Rn_inv, rho, L, Q, gamma0)
% PRD_design - 算法3：PRD 波形设计总体流程（论文 Algorithm 3）
%
% 步骤：
%   1. 预计算系数矩阵
%   2. 求解初始化问题 P0，判断高/低 SINR
%   3. 调用算法1（高SINR）或算法2（低SINR）优化 Rs
%   4. 从 Rs 合成波形 s*
%
% 输入：
%   t0, Rt  - TIR 均值和协方差
%   Rn_inv  - 干扰协方差逆矩阵（(L+Q-1) x (L+Q-1)）
%   rho     - PAR 约束上限
%   L, Q    - 波形长度和目标单元数
%   gamma0  - SINR 阈值（线性值）
%
% 输出：
%   s_star    - 合成的最优波形（L x 1）
%   Rs_opt    - 优化的波形协方差矩阵（L x L）
%   f_Rs      - f(Rs)，矩阵层面的 PRD 值
%   f_s       - f(s*)，实际波形的 PRD 值
%   f_vals    - 每次迭代的 Dinkelbach 辅助目标值数组
%   sinr_case - 'high' 或 'low'

fprintf('预计算系数矩阵...');
coeffs = precompute_coeffs(t0, Rt, Rn_inv, Q, L);
fprintf('完成\n');

% Step 1：求解初始化问题 P0
fprintf('求解初始化问题 P0...');
[Rs0, f0, obj0] = solve_init_SDP(t0, Rt, rho, L, Q, gamma0, coeffs);
fprintf('完成，obj0 = %.4f\n', obj0);

% Step 2：判断 SINR 情况（论文 Section III.C）
if obj0 >= 0
    sinr_case = 'high';
    fprintf('判断为高SINR情况，调用算法1（Dinkelbach+SDP）...\n');
    [Rs_opt, f_vals] = solve_high_SINR(t0, Rt, rho, L, Q, gamma0, coeffs, Rs0, f0);
else
    sinr_case = 'low';
    fprintf('判断为低SINR情况，调用算法2（MM+Dinkelbach+SDP）...\n');
    [Rs_opt, f_vals] = solve_low_SINR(t0, Rt, rho, L, Q, gamma0, coeffs, Rs0, f0);
end
fprintf('优化完成，共 %d 次迭代\n', length(f_vals));

% Step 3：计算 f(Rs)
M_opt = compute_M(Rs_opt, coeffs.F, Q);
[~, ~, f_Rs] = compute_EV_from_M(M_opt, t0, Rt, gamma0);

% Step 4：波形合成
fprintf('波形合成（随机化）...');
[s_star, f_s] = synthesize_waveform(Rs_opt, t0, Rt, coeffs, Q, L, gamma0, rho);
fprintf('完成，f(s*) = %.4f\n', f_s);
end
