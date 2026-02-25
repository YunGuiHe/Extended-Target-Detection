function p = get_target_params(target_id)
% get_target_params - 返回论文中两个目标的 TIR 参数
%
% 输入：
%   target_id - 1（战斗机）或 2（客机）
% 输出：
%   p.L   - 波形采样点数（60）
%   p.Q   - 目标距离单元数
%   p.t0  - TIR 均值向量（Q x 1）
%   p.Rt  - TIR 协方差矩阵（Q x Q）

p.L = 60;

if target_id == 1
    % 目标1：战斗机，15m，Q = 10 个距离单元
    p.Q = 10;
    q = (1:p.Q)';
    p.t0 = (cos(0.6*pi*q - 1.8*pi) + 1) .* exp(-abs(q - 2) + 1i*pi*(q - 4)/6);
    p.Rt = zeros(p.Q, p.Q);
    for m = 1:p.Q
        for n = 1:p.Q
            p.Rt(m,n) = 1.2*exp(-(m-n)^2) + 0.2*exp(1i*2*(m-n));
        end
    end

elseif target_id == 2
    % 目标2：客机，38m，Q = 25 个距离单元
    p.Q = 25;
    p.t0 = 0.1 * ones(p.Q, 1);
    p.Rt = zeros(p.Q, p.Q);
    for m = 1:p.Q
        for n = 1:p.Q
            p.Rt(m,n) = 0.2*exp(1i*2*(m-n)) + 0.1*exp(1i*4*(m-n)) + 0.3*exp(1i*6*(m-n));
        end
    end

else
    error('target_id 必须为 1 或 2');
end
end
