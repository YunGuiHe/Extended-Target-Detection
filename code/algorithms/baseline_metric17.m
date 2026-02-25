function [s_ref, Rs_ref, f_Rs_ref, f_s_ref] = baseline_metric17(t0, Rt, Rn_inv, rho, L, Q, gamma0)
% baseline_metric17 - 文献[17]对比基准：只最大化 E[gamma]（不考虑方差）
%
% 优化问题（单次线性 SDP，无迭代）：
%   max  tr(Rt*M) + t0^H*M*t0  [= E[gamma]]
%   s.t. tr(Rs) = 1, Rs >= 0, Rs(i,i) <= rho
%
% 波形合成时按 E[gamma]（而非 PRD 指标 f）选择最优波形
%
% 输出：
%   s_ref      - 文献[17]方法合成的波形（L x 1）
%   Rs_ref     - 波形协方差矩阵
%   f_Rs_ref   - 该波形对应的 PRD 指标 f（用于对比）
%   f_s_ref    - 合成波形的 PRD 指标

coeffs = precompute_coeffs(t0, Rt, Rn_inv, Q, L);
F      = coeffs.F;

% 构造 E[gamma] 线性系数矩阵
Z_E = zeros(L, L);
for m = 1:Q
    for n = 1:Q
        coeff = Rt(n,m) + conj(t0(m)) * t0(n);
        Z_E   = Z_E + coeff * F{m,n};
    end
end

cvx_begin sdp quiet
    variable Rs_var(L, L) hermitian semidefinite
    maximize(real(trace(Z_E * Rs_var)))
    subject to
        trace(Rs_var) == 1
        real(diag(Rs_var)) <= rho * ones(L, 1)
cvx_end

Rs_ref = (Rs_var + Rs_var') / 2;

% 计算 PRD 指标（用于和 PRD 方法的指标作对比）
M_ref = compute_M(Rs_ref, F, Q);
[~, ~, f_Rs_ref] = compute_EV_from_M(M_ref, t0, Rt, gamma0);

% 合成波形：按 E[gamma] 选择最优波形（而非 PRD）
[s_ref, f_s_ref] = synthesize_waveform_by_Egamma(Rs_ref, t0, Rt, coeffs, Q, L, gamma0, rho);
end
