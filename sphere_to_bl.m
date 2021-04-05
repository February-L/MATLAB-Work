function [r_bl, theta_bl] = sphere_to_bl(r_sphere, theta_sphere, a)
% 将球坐标转换为Kerr时空转动参数为a的B-L坐标（两者不同的仅为r、theta坐标的定义）。
% 使用自然单位制：c = G = 1。
%
% function [r_bl, theta_bl] = sphere_to_bl(r_sphere, theta_sphere, a)
%     r_sphere、theta_sphere分别指定球坐标的r、theta值，应为具有相同长度的数组；
%     a指定B-L坐标相关的Kerr时空转动参数a的值，a与r同量纲。
%     返回给定球坐标对应的B-L坐标r、theta值，r_bl与theta_bl均为与r_sphere、theta_sphere同长度的数组。
%     B-L坐标r = 0处的theta_bl返回值的在区间[0, pi/2]内。
%     只保证对非负实值的r_sphere、实值的theta_sphere以及非零实值的a给出正确的返回结果。
    r_bl = sqrt(0.5 .* ( ...
               r_sphere.^2.0 - a.^2.0 + ...
               sqrt((r_sphere.^2.0 - a.^2.0).^2.0 + 4.0 .* (a .* r_sphere .* cos(theta_sphere)).^2.0) ...
           ));
    theta_bl = acos(r_sphere .* cos(theta_sphere) ./ r_bl);
    sp_indexes = find(isnan(theta_bl) | isinf(theta_bl));
    theta_bl(sp_indexes) = asin(r_sphere(sp_indexes) ./ abs(a));
end
