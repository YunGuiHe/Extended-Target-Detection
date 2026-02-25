function M = compute_M(Rs, F, Q)
% compute_M - 数值计算 M 矩阵（用于 CVX 外部评估）
%
% M(m,n) = trace(F{m,n} * Rs)
%
% 输入：
%   Rs - 波形协方差矩阵（L x L）
%   F  - 预计算系数矩阵 cell{m,n}（每个为 L x L）
%   Q  - 目标距离单元数
% 输出：
%   M  - Q x Q 矩阵

M = zeros(Q, Q);
for m = 1:Q
    for n = 1:Q
        M(m, n) = trace(F{m,n} * Rs);
    end
end
end
