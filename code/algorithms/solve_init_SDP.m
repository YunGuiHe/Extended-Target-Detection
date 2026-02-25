function [Rs_opt, f_opt, obj_val] = solve_init_SDP(t0, Rt, rho, L, Q, gamma0, coeffs)
Z_num = coeffs.Z_num;
F = coeffs.F;
cvx_begin sdp quiet
    variable Rs_init(L, L) hermitian semidefinite
    maximize( real(trace(Z_num * Rs_init)) - gamma0 )
    subject to
        trace(Rs_init) == 1;
        real(diag(Rs_init)) <= rho * ones(L, 1);
cvx_end
obj_val = cvx_optval;
Rs_opt = (Rs_init + Rs_init') / 2;
M = compute_M(Rs_opt, F, Q);
[~, ~, f_opt] = compute_EV_from_M(M, t0, Rt, gamma0);
end
