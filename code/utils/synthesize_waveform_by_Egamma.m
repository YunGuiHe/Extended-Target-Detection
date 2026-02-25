function [s_star, f_star] = synthesize_waveform_by_Egamma(Rs, t0, Rt, coeffs, Q, L, gamma0, rho, N_rand)
% synthesize_waveform_by_Egamma - 从 Rs 合成波形，按 E[gamma] 选择最优（用于[17]基准方法）
%
% CMC 情况下使用相位提取法，选择准则为 E[gamma] = tr(Rt*M) + t0'*M*t0

if nargin < 9, N_rand = 1000; end
F = coeffs.F;

is_cmc = (abs(rho - 1/L) < 1e-10);

% 特征分解
[V_eig, D_eig] = eig(Rs);
[d_sorted, idx_sorted] = sort(real(diag(D_eig)), 'descend');
d_sorted = max(d_sorted, 0);
V_sorted = V_eig(:, idx_sorted);

is_rank1 = (sum(d_sorted) > 1e-12 && d_sorted(1) / sum(d_sorted) > 0.9999);

if is_rank1 && ~is_cmc
    s_star = V_sorted(:, 1);
    s_star = s_star / norm(s_star);
    s_star = par_clip(s_star, rho);
    M = compute_M(s_star * s_star', F, Q);
    [~,~,f_star] = compute_EV_from_M(M, t0, Rt, gamma0);
    return;
end

%% 多策略合成（按 E[gamma] 选择）
E_best = -Inf;
s_star = zeros(L, 1);

% 策略1：主特征向量相位
v1 = V_sorted(:, 1);
s_c = exp(1i * angle(v1)) / sqrt(L);
if ~is_cmc, s_c = par_clip(s_c, rho); end
s_c = s_c / norm(s_c);
[~, s_star, E_best] = eval_E_and_update(s_c, F, Q, t0, Rt, s_star, E_best);

% 策略2：前几个特征向量
n_top = min(5, sum(d_sorted > 1e-6));
for k = 1:n_top
    vk = V_sorted(:, k);
    s_c = exp(1i * angle(vk)) / sqrt(L);
    if ~is_cmc, s_c = par_clip(s_c, rho); end
    s_c = s_c / norm(s_c);
    [~, s_star, E_best] = eval_E_and_update(s_c, F, Q, t0, Rt, s_star, E_best);
end

% 策略3：加权叠加
for trial = 1:min(5, n_top)
    weights = d_sorted(1:n_top) .* (randn(n_top, 1) + 1i*randn(n_top, 1));
    v_combo = V_sorted(:, 1:n_top) * weights;
    s_c = exp(1i * angle(v_combo)) / sqrt(L);
    if ~is_cmc, s_c = par_clip(s_c, rho); end
    s_c = s_c / norm(s_c);
    [~, s_star, E_best] = eval_E_and_update(s_c, F, Q, t0, Rt, s_star, E_best);
end

% 策略4：随机化
Rs_sym = (Rs + Rs') / 2 + 1e-9 * eye(L);
try
    Rs12 = chol(Rs_sym, 'lower');
catch
    [V_tmp, D_tmp] = eig(Rs_sym);
    d_tmp = max(real(diag(D_tmp)), 0);
    Rs12 = V_tmp * diag(sqrt(d_tmp));
end

for i = 1:N_rand
    v = (randn(L, 1) + 1i * randn(L, 1)) / sqrt(2);
    s_c = Rs12 * v;
    if norm(s_c) < 1e-12, continue; end
    
    if is_cmc
        s_c = exp(1i * angle(s_c)) / sqrt(L);
    else
        s_c = s_c / norm(s_c);
        s_c = par_clip(s_c, rho);
    end
    s_c = s_c / norm(s_c);
    
    [~, s_star, E_best] = eval_E_and_update(s_c, F, Q, t0, Rt, s_star, E_best);
end

% 计算最终的 PRD 指标
M_star = compute_M(s_star * s_star', F, Q);
[~,~,f_star] = compute_EV_from_M(M_star, t0, Rt, gamma0);
end

function s = par_clip(s, rho)
    for iter = 1:5
        over = abs(s).^2 > rho + 1e-9;
        if ~any(over), break; end
        s(over) = s(over) .* sqrt(rho) ./ abs(s(over));
        s = s / norm(s);
    end
end

function [E_c, s_best, E_best] = eval_E_and_update(s_c, F, Q, t0, Rt, s_best, E_best)
    M_c = compute_M(s_c * s_c', F, Q);
    E_c = real(trace(Rt * M_c)) + real(t0' * M_c * t0);
    if E_c > E_best
        E_best = E_c;
        s_best = s_c;
    end
end
