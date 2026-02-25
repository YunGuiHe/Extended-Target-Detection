function [Rs_opt, f_vals] = solve_high_SINR(t0, Rt, rho, L, Q, gamma0, coeffs, Rs0, f0)
% solve_high_SINR - 算法1：Dinkelbach + SDP 的高 SINR 波形设计
%
% 将 f(Rs) = (E[gamma]-gamma0)/D[gamma] 通过 Dinkelbach 迭代转化为：
%   min  -(E[gamma]-gamma0) + f^(k) * D[gamma]
% 其中 D[gamma] = sqrt(V[gamma])，通过 SOC 约束处理
%
% 方差展开为 Frobenius 范数平方：
%   V = ||Rt^{1/2} * M * K||_F^2
%   其中 K = [Rt^{1/2}, sqrt(2)*t0]
%   元素级别：z_{ij} = tr(Z_B{i,j} * Rs)
%   V = sum |z_{ij}|^2 = norm([Re(z); Im(z)])^2

Z_num = coeffs.Z_num;
Z_B = coeffs.Z_B;
tol = 1e-4;
max_iter = 50;
N_b = Q * (Q + 1);

% 预计算 Z_B 的厄米部分和反厄米部分，用于分离实虚部
Z_B_h = cell(Q, Q+1);
Z_B_s = cell(Q, Q+1);
for i = 1:Q
    for j = 1:Q+1
        Z_B_h{i,j} = (Z_B{i,j} + Z_B{i,j}') / 2;
        Z_B_s{i,j} = (Z_B{i,j} - Z_B{i,j}') / (2i);
    end
end

Rs_cur = Rs0;
f_cur = f0;
f_vals = [];

for k = 1:max_iter
    f_prev = f_cur;
    fp = max(f_prev, 1e-12);

    cvx_begin sdp quiet
        variable Rs_var(L, L) hermitian semidefinite
        expression b_re(N_b, 1)
        expression b_im(N_b, 1)

        idx = 0;
        for ii = 1:Q
            for jj = 1:Q+1
                idx = idx + 1;
                b_re(idx) = trace(Z_B_h{ii,jj} * Rs_var);
                b_im(idx) = trace(Z_B_s{ii,jj} * Rs_var);
            end
        end

        % E[gamma] - gamma0
        num_expr = real(trace(Z_num * Rs_var)) - gamma0;

        % min -(E-gamma0) + f^(k) * sqrt(V)
        % sqrt(V) = norm([b_re; b_im])
        minimize(-num_expr + fp * norm([b_re; b_im]))

        subject to
            trace(Rs_var) == 1;
            real(diag(Rs_var)) <= rho * ones(L, 1);
            num_expr >= 0;
    cvx_end

    if ~strcmp(cvx_status, 'Solved') && ~strcmp(cvx_status, 'Inaccurate/Solved')
        warning('CVX failed at iter %d: %s', k, cvx_status);
        break;
    end

    Rs_new = (Rs_var + Rs_var') / 2;
    M_new = compute_M(Rs_new, coeffs.F, Q);
    [E_new, V_new, f_new] = compute_EV_from_M(M_new, t0, Rt, gamma0);

    % Dinkelbach 辅助目标值（收敛到0）
    aux_val = (E_new - gamma0) - f_prev * sqrt(max(V_new, 0));
    f_vals(end+1) = aux_val;  %#ok<AGROW>

    Rs_cur = Rs_new;
    f_cur = f_new;

    if abs(aux_val) < tol
        break;
    end
    if k > 1 && abs(f_vals(end) - f_vals(end-1)) / (abs(f_vals(end-1)) + 1e-12) < tol
        break;
    end
end

Rs_opt = Rs_cur;
end
