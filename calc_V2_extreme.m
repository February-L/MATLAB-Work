function extreme_r = calc_V2_extreme(M, g, L, is_timelike)
    % function extreme_r = calc_V2_extreme(M, g, L, is_timelike)
    %     根据Bardeen时空的引力源参数（质量M与正则化参数g）与物体关于转动角phi的守恒量L计算有效势极值点位置。
    %     使用自然单位制：c = G = 1, mu_0 = 4*pi。
    %     is_timelike为指定物体运动是否类时（否则视为类光）的布尔值参数。
    %     当is_timelike为假时，类光运动的特点使得L的值是无所谓的。
    %     返回由极值点r坐标（由内到外）组成的行数组（元素数为极值点数）。
    if (is_timelike)
        tmp_a = g^2.0; tmp_b = M^2.0 / L^4.0; tmp_c = 2.0 * g^2.0 - 3.0 * L^2.0;
        derivative_roots = sqrt(roots([-tmp_b;
                                       1.0 + 2.0 * tmp_b * tmp_c;
                                       5.0 * tmp_a - tmp_b * tmp_c^2.0;
                                       10.0 * tmp_a^2.0;
                                       10.0 * tmp_a^3.0;
                                       5.0 * tmp_a^4.0;
                                       tmp_a^5.0]));
    else
        tmp_a = g^2.0; tmp_b = M^2.0;
        derivative_roots = sqrt(roots([1.0;
                                       5.0 * tmp_a - 9.0 * tmp_b;
                                       10.0 * tmp_a^2.0;
                                       10.0 * tmp_a^3.0;
                                       5.0 * tmp_a^4.0;
                                       tmp_a^5.0]));
    end
    root_count = length(derivative_roots);
    extreme_r = zeros(1, 3);
    for i = 1: 1: root_count
        if ((imag(derivative_roots(i)) == 0.0) && (real(derivative_roots(i)) > 0.0))
            if (extreme_r(1) == 0.0)
                extreme_r(1) = derivative_roots(i);
            elseif (extreme_r(2) == 0.0)
                extreme_r(2) = derivative_roots(i);
            else
                extreme_r(3) = derivative_roots(i);
            end
        end
    end
    extreme_r = unique(extreme_r(extreme_r ~= 0.0));
    extreme_r = sort(extreme_r, 'ascend');
end
