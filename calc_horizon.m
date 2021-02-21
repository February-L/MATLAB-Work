function horizon_r = calc_horizon(M, g)
    % function horizon_r = calc_horizon(M, g)
    %     根据Bardeen时空的引力源参数（质量M与正则化参数g）计算视界半径。
    %     使用自然单位制：c = G = 1, mu_0 = 4*pi。
    %     返回由视界半径（由内到外）组成的行数组（元素数为视界数），无视界即返回空数组。
    %     (g / M)^2 < 16/27时才有视界，且g不等于0时为双视界(g = 0则退化至Schwarzschild时空，r = 0处形成奇点)。
    metric_roots = sqrt(roots([1.0; 3.0 * g^2.0 - 4.0 * M^2.0; 3.0 * g^4.0; g^6.0]));
    horizon_r = zeros(1, 2);
    for i = 1: 1: 3
        if ((imag(metric_roots(i)) == 0.0) && (real(metric_roots(i)) > 0.0))
            if (horizon_r(1) == 0.0)
                horizon_r(1) = metric_roots(i);
            else
                horizon_r(2) = metric_roots(i);
            end
        end
    end
    horizon_r = unique(horizon_r(horizon_r ~= 0.0));
    horizon_r = sort(horizon_r, 'ascend');
end
