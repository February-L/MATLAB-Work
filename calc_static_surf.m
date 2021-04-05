function r_static = calc_static_surf(the_theta, M, a)
% 计算Kerr时空静止界限在给定方位（以B-L坐标theta表示）的位置（以B-L坐标r表示）。
% 使用自然单位制：c = G = 1。
%
% function r_static = calc_static_surf(the_theta, M, a)
%     the_theta指定方位，为B-L坐标theta的值（数组）。
%     M、a分别指定引力源的质量与转动参数值，a与M同量纲。
%     返回由两静止界限在给定方位（序列）的B-L坐标r值构成的2行矩阵r_static，靠内静界在第1行，若该方位上不存在静止界限则对应值为NaN。
    if (M <= 0.0)
        r_static = zeros(2, numel(the_theta)) + NaN;
        return;
    end
    the_theta = (the_theta(:))';
    r_static = [M - sqrt(M^2.0 - (a .* cos(the_theta)).^2.0);
                M + sqrt(M^2.0 - (a .* cos(the_theta)).^2.0)];
    r_static(imag(r_static) ~= 0.0) = NaN;
end
