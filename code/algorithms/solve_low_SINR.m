function [Rs_opt, f_vals] = solve_low_SINR(t0, Rt, rho, L, Q, gamma0, coeffs, Rs0, f0)
% solve_low_SINR - 算法2：MM（内层）+ Dinkelbach（外层）的低 SINR 波形设计
%
% 适用条件：所有可行波形均有 E[gamma] < gamma0（初始化值 < 0）
%
% 外层 Dinkelbach 迭代，内层 MM 迭代：
%   内层 MM 子问题（每步为线性 SDP，论文公式 18）：
%     max  tr[(Rt + t0*t0^H - f^(k-1)*G^(k,v)) * M]
%     s.t. tr(Rs)=1, Rs>=0, Rs(i,i)<=rho
%
% 梯度矩阵（论文公式 17）：
%   G^(k) = Rt*M^(k)*(Rt + t0*t0^H) / sqrt(tr(Rt*M^(k)*Rt*M^(k)) + t0^H*M^(k)*Rt*M^(k)*t0)
%   注意：分母中第二项系数为1（非2），这与方差公式不同
%
% 输出：
%   Rs_opt - 次优的波形协方差矩阵
%   f_vals - 外层每次迭代后的 Dinkelbach 辅助目标值

F  = coeffs.F;

tol_outer  = 1e-4;
tol_inner  = 1e-5;
max_outer  = 30;
max_inner  = 30;

Rs_cur = Rs0;
f_cur  = f0;
f_vals = [];

for k = 1:max_outer
    f_prev  = f_cur;
    Rs_mm   = Rs_cur;  % 内层 MM 的初始点

    for v = 1:max_inner
        % 计算当前 M
        M_cur = compute_M(Rs_mm, F, Q);
        RtM   = Rt * M_cur;

        % 计算真正的方差 V
        term1 = real(trace(RtM * RtM));
        term2 = real(t0' * M_cur * RtM * t0);
        D_G = sqrt(max(term1 + 2 * term2, 1e-12));

        % 梯度矩阵 G = (Rt*M*Rt + Rt*M*t0*t0' + t0*t0'*M*Rt) / D_G (严格的解析梯度)
        Rt_plus = Rt + t0 * t0';
        G_numerator = RtM * Rt + RtM * (t0 * t0') + (t0 * t0') * M_cur * Rt;
        G = G_numerator / D_G;

        % 构造内层线性 SDP 系数矩阵 Z_C
        % obj = tr[(Rt + t0*t0^H - f_prev*G) * M] = tr(Z_C * Rs)
        C_coeff = Rt_plus - f_prev * G;   % Q x Q

        Z_C = zeros(L, L);
        for mm = 1:Q
            for nn = 1:Q
                % tr(C*M) = sum_{m,n} C(m,n)*M(n,m) = sum C(n,m)*tr(F{m,n}*Rs)
                Z_C = Z_C + C_coeff(nn, mm) * F{mm, nn};
            end
        end

        % 求解内层线性 SDP
        cvx_begin sdp quiet
            variable Rs_var(L, L) hermitian semidefinite
            maximize(real(trace(Z_C * Rs_var)))
            subject to
                trace(Rs_var) == 1
                real(diag(Rs_var)) <= rho * ones(L, 1)
        cvx_end

        if ~strcmp(cvx_status, 'Solved') && ~strcmp(cvx_status, 'Inaccurate/Solved')
            warning('CVX 内层求解失败（外层%d，内层%d）：%s', k, v, cvx_status);
            break;
        end

        Rs_new = (Rs_var + Rs_var') / 2;
        M_new  = compute_M(Rs_new, F, Q);

        % 内层收敛判断（目标函数值变化）
        obj_old = real(trace(C_coeff * M_cur));
        obj_new = real(trace(C_coeff * M_new));
        Rs_mm   = Rs_new;

        if abs(obj_new - obj_old) / (abs(obj_old) + 1e-12) < tol_inner
            break;
        end
    end

    % 外层：更新 f 值
    Rs_cur = Rs_mm;
    M_cur  = compute_M(Rs_cur, F, Q);
    [E_new, V_new, f_new] = compute_EV_from_M(M_cur, t0, Rt, gamma0);
    f_cur  = f_new;

    % Dinkelbach 辅助目标值
    if V_new > 0
        aux_val = (E_new - gamma0) - f_prev * sqrt(V_new);
    else
        aux_val = E_new - gamma0;
    end
    f_vals(end+1) = aux_val;   %#ok<AGROW>

    % 外层收敛判断
    if k > 1 && abs(f_vals(end) - f_vals(end-1)) / (abs(f_vals(end-1)) + 1e-12) < tol_outer
        break;
    end
end

Rs_opt = Rs_cur;
end
