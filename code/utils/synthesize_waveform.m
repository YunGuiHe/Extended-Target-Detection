function [s_star, f_star] = synthesize_waveform(Rs, t0, Rt, coeffs, Q, L, gamma0, rho, N_rand)
% synthesize_waveform - 从波形协方差矩阵 Rs 合成满足 PAR 约束的波形 s*
%
% CMC 情况(rho ≈ 1/L)下使用相位提取法：
%   1. 对 Rs 特征分解，提取主特征向量的相位
%   2. 随机化时直接投影到恒模流形（保留相位，固定幅度）
%
% 一般 PAR 情况：
%   1. 若 Rs 近似秩-1，直接做特征分解
%   2. 否则：随机化 N_rand 个候选波形，选择 f(si) 最大的

if nargin < 9, N_rand = 1000; end
F = coeffs.F;

is_cmc = (abs(rho - 1/L) < 1e-10);  % 判断是否为恒模约束

% 特征分解
[V_eig, D_eig] = eig(Rs);
[d_sorted, idx_sorted] = sort(real(diag(D_eig)), 'descend');
d_sorted = max(d_sorted, 0);
V_sorted = V_eig(:, idx_sorted);

% 检查 Rs 是否接近秩-1
is_rank1 = (sum(d_sorted) > 1e-12 && d_sorted(1) / sum(d_sorted) > 0.9999);

if is_rank1 && ~is_cmc
    % 近似秩-1 且非CMC：取主特征向量
    s_star = V_sorted(:, 1);
    s_star = s_star / norm(s_star);
    s_star = par_clip(s_star, rho);
    M = compute_M(s_star * s_star', F, Q);
    [~,~,f_star] = compute_EV_from_M(M, t0, Rt, gamma0);
    return;
end

%% CMC 或非秩-1情况：多策略合成
f_best = -Inf;
s_star = zeros(L, 1);

% 策略1：主特征向量相位提取（CMC最重要的策略）
v1 = V_sorted(:, 1);
s_c = exp(1i * angle(v1)) / sqrt(L);  % 提取相位，恒模
if ~is_cmc
    s_c = par_clip(s_c, rho);
end
s_c = s_c / norm(s_c);
[f_c, s_star, f_best] = eval_and_update(s_c, F, Q, t0, Rt, gamma0, s_star, f_best);

% 策略2：前几个特征向量的加权组合相位提取
n_top = min(5, sum(d_sorted > 1e-6));
for k = 1:n_top
    vk = V_sorted(:, k);
    s_c = exp(1i * angle(vk)) / sqrt(L);
    if ~is_cmc
        s_c = par_clip(s_c, rho);
    end
    s_c = s_c / norm(s_c);
    [~, s_star, f_best] = eval_and_update(s_c, F, Q, t0, Rt, gamma0, s_star, f_best);
end

% 策略3：加权叠加后提取相位
for trial = 1:min(5, n_top)
    weights = d_sorted(1:n_top) .* (randn(n_top, 1) + 1i*randn(n_top, 1));
    v_combo = V_sorted(:, 1:n_top) * weights;
    s_c = exp(1i * angle(v_combo)) / sqrt(L);
    if ~is_cmc
        s_c = par_clip(s_c, rho);
    end
    s_c = s_c / norm(s_c);
    [~, s_star, f_best] = eval_and_update(s_c, F, Q, t0, Rt, gamma0, s_star, f_best);
end

% 策略4：随机化（CMC投影 vs PAR裁剪）
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
        % CMC投影：保留相位，设置恒模幅度
        s_c = exp(1i * angle(s_c)) / sqrt(L);
    else
        s_c = s_c / norm(s_c);
        s_c = par_clip(s_c, rho);
    end
    s_c = s_c / norm(s_c);
    
    [~, s_star, f_best] = eval_and_update(s_c, F, Q, t0, Rt, gamma0, s_star, f_best);
end

f_star = f_best;
end

%% 辅助函数
function s = par_clip(s, rho)
    for iter = 1:5
        over = abs(s).^2 > rho + 1e-9;
        if ~any(over), break; end
        s(over) = s(over) .* sqrt(rho) ./ abs(s(over));
        s = s / norm(s);
    end
end

function [f_c, s_best, f_best] = eval_and_update(s_c, F, Q, t0, Rt, gamma0, s_best, f_best)
    M_c = compute_M(s_c * s_c', F, Q);
    [~, ~, f_c] = compute_EV_from_M(M_c, t0, Rt, gamma0);
    if f_c > f_best
        f_best = f_c;
        s_best = s_c;
    end
end
