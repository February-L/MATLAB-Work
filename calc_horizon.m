function r_horizon = calc_horizon(M, a)
% 计算Kerr时空视界的位置（以B-L坐标r表示）。
% 使用自然单位制：c = G = 1。
%
% function r_horizon = calc_horizon(M, a)
%     M、a分别指定引力源的质量与转动参数值，a与M同量纲。
%     返回由两视界的B-L坐标r值构成的2元素数组r_horizon，靠内视界在前，若视界不存在则对应值为NaN。
%     在|a| = M而使两视界于r = M处完全重合时，返回的r_horizon将为[M; NaN]。
    a = abs(a);
    if (a < M)
        r_horizon = M + sqrt(M^2.0 - a^2.0) .* [-1.0; 1.0];
    elseif ((a == M) && (M ~= 0.0))
        r_horizon = [M; NaN];
    else
        r_horizon = [NaN; NaN];
    end
end
