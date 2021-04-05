function effp_val = eff_potential_r(r, is_upper, is_timelike, L, M, a, clean_imag_err)
% 计算Kerr时空外部赤道面（B-L坐标theta = pi/2）上类时或类光运动的上（下）有效势V+（V-）在给定位置的值。
% 后缀r代表real，复值返回为NaN。
% 使用自然单位制：c = G = 1。
%
% function effp_val = eff_potential_r(r, is_upper, is_timelike, L, M, a, clean_imag_err)
%     r指定位置，为B-L坐标r的值；
%     is_upper指定所求是否为上有效势，否则为下有效势；
%     is_timelike指定所求是否为类时有效势，否则为类光有效势；
%     L指定运动关于B-L坐标phi的守恒量的值；
%     M、a分别指定引力源的质量与转动参数值，a与M同量纲；
%     clean_imag_err指定是否将（视界附近）纯由误差导致的虚部抹去。
%     返回所指定的有效势在指定位置（数组）的值（数组）effp_val，复值返回为NaN。
    upper_val = double(logical(is_upper)) * 2.0 - 1.0;
    timelike_val = double(logical(is_timelike));
    M_dbl = 2.0 * M; a_sqr = a ^ 2.0;
    effp_val = 1.0 ./ (r.^2.0 + a_sqr + M_dbl .* a_sqr ./ r) .* ( ...
                   M_dbl .* a .* L ./ r + upper_val .* sqrt( ...
                       (r.^2.0 - M_dbl .* r + a_sqr) .* (L^2.0 + timelike_val .* ( ...
                           r.^2.0 + a_sqr + M_dbl .* a_sqr ./ r ...
                       )) ...
                   ) ...
               );
    if (a ~= 0.0)
        effp_val(r == 0.0) = L / a;
     end
    if (clean_imag_err && (abs(a) <= M))
        r_horizon = M + sqrt(M^2.0 - a^2.0) .* [-1.0, 1.0];
        outer_imag_i = find(((r <= r_horizon(1)) | (r >= r_horizon(2))) & (imag(effp_val) ~= 0.0));
        effp_val(outer_imag_i) = real(effp_val(outer_imag_i));
    end
    effp_val(imag(effp_val) ~= 0.0) = NaN;
end
