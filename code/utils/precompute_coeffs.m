function coeffs = precompute_coeffs(t0, Rt, Rn_inv, Q, L)
% precompute_coeffs - 预计算优化所需的各种系数矩阵
%
% 计算：
%   F{m,n}  = J_m' * Rn_inv * J_n   （M矩阵的线性表达系数）
%   Z_num   = sum (Rt(m,n)+t0(m)'*t0(n)) * F{m,n}  （E[gamma]的线性系数）
%   Z_B{i,j} = sum_m Rt12(i,m) * G{m,j}  （方差的Frobenius范数展开系数）
%   其中 G{m,j} = sum_n ColK(n,j) * F{m,n}
%        ColK = [Rt12, sqrt(2)*t0]
%
% 方差展开：V[gamma] = ||Rt12 * M * ColK||_F^2

% 使用特征分解计算 Rt^{1/2}，确保半正定
[V_eig, D_eig] = eig(Rt);
d = real(diag(D_eig));
d = max(d, 0);   % 截断负特征值
Rt12 = V_eig * diag(sqrt(d)) * V_eig';

ColK = [Rt12, sqrt(2)*t0];

J_cells = cell(Q, 1);
for q = 0:Q-1
    J_cells{q+1} = compute_J(q, L, Q);
end

RnJ = cell(Q, 1);
for n = 1:Q
    RnJ{n} = Rn_inv * J_cells{n};
end

F = cell(Q, Q);
for m = 1:Q
    for n = 1:Q
        F{m,n} = J_cells{m}' * RnJ{n};
    end
end

Z_num = zeros(L, L);
for m = 1:Q
    for n = 1:Q
        Z_num = Z_num + (Rt(n,m) + conj(t0(m))*t0(n)) * F{m,n};
    end
end

G = cell(Q, Q+1);
for m = 1:Q
    for j = 1:Q+1
        Gmj = zeros(L, L);
        for n = 1:Q
            Gmj = Gmj + ColK(n,j) * F{m,n};
        end
        G{m,j} = Gmj;
    end
end

Z_B = cell(Q, Q+1);
for ii = 1:Q
    for j = 1:Q+1
        Zij = zeros(L, L);
        for m = 1:Q
            Zij = Zij + Rt12(ii,m) * G{m,j};
        end
        Z_B{ii,j} = Zij;
    end
end

coeffs.J_cells = J_cells;
coeffs.F = F;
coeffs.Z_num = Z_num;
coeffs.Z_B = Z_B;
coeffs.Rt12 = Rt12;
coeffs.ColK = ColK;
end
