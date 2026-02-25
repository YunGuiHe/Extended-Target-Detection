clear; clc;
addpath('utils'); addpath('algorithms');
p = get_target_params(1);
Q = p.Q; L = p.L; t0 = p.t0; Rt = p.Rt;
sigma2 = 10^(-10/10);
Rn_inv = (1/sigma2) * eye(L + Q - 1);
coeffs = precompute_coeffs(t0, Rt, Rn_inv, Q, L);
F = coeffs.F;

% Test Random M
s = (randn(L,1) + 1i*randn(L,1))/sqrt(2);
s = s / norm(s);
M = compute_M(s*s', F, Q);
[E_true, ~, ~] = compute_EV_from_M(M, t0, Rt, 1);

% Test Z_num with Rt(m,n)
Z_wrong = zeros(L,L);
Z_right = zeros(L,L);
for m=1:Q
    for n=1:Q
        Z_wrong = Z_wrong + (Rt(m,n) + conj(t0(m))*t0(n)) * F{m,n};
        Z_right = Z_right + (Rt(n,m) + conj(t0(m))*t0(n)) * F{m,n};
    end
end

E_wrong = real(trace(Z_wrong * (s*s')));
E_right = real(trace(Z_right * (s*s')));

fprintf('True E_val = %f\n', E_true);
fprintf('E_val from Rt(m,n) = %f\n', E_wrong);
fprintf('E_val from Rt(n,m) = %f\n', E_right);
