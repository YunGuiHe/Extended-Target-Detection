function J = compute_J(q, L, Q)
% compute_J - 构造移位矩阵 J_q，大小 (L+Q-1) x L
% J_q(m,n) = 1 当 m-n = q（1-indexed），否则为 0
%
% 输入：
%   q  - 移位量（0-indexed，范围 0 到 Q-1）
%   L  - 波形采样点数
%   Q  - 目标占据的距离单元数
% 输出：
%   J  - (L+Q-1) x L 移位矩阵

N = L + Q - 1;
J = zeros(N, L);
for n = 1:L
    m = n + q;
    if m >= 1 && m <= N
        J(m, n) = 1;
    end
end
end
