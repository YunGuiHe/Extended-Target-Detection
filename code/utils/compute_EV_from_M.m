function [E_val, V_val, f_val] = compute_EV_from_M(M, t0, Rt, gamma0)
% compute_EV_from_M - 由 M 矩阵计算 E[gamma]、V[gamma] 和 PRD 指标 f
%
% 论文公式（Proposition 1，Appendix A）：
%   E[gamma] = tr(Rt*M) + t0^H*M*t0
%   V[gamma] = tr(Rt*M*Rt*M) + 2*t0^H*M*Rt*M*t0  （第二项系数为 2）
%   f(s)     = (E[gamma] - gamma0) / sqrt(V[gamma])
%
% 输入：
%   M      - Q x Q 矩阵
%   t0     - TIR 均值（Q x 1）
%   Rt     - TIR 协方差（Q x Q）
%   gamma0 - SINR 阈值（线性值，非 dB）
% 输出：
%   E_val, V_val, f_val

RtM = Rt * M;
E_val = real(trace(RtM)) + real(t0' * M * t0);
V_val = real(trace(RtM * RtM)) + 2 * real(t0' * M * RtM * t0);

if V_val <= 1e-12
    f_val = sign(E_val - gamma0) * Inf;
else
    f_val = (E_val - gamma0) / sqrt(V_val);
end
end
