function [r_sphere, theta_sphere] = bl_to_sphere(r_bl, theta_bl, a)
% 将Kerr时空转动参数为a的B-L坐标转换为球坐标（两者不同的仅为r、theta坐标的定义）。
% 使用自然单位制：c = G = 1。
%
% function [r_sphere, theta_sphere] = bl_to_sphere(r_bl, theta_bl, a)
%     r_bl、theta_bl分别指定B-L坐标的r、theta值，应为具有相同长度的数组；
%     a指定B-L坐标相关的Kerr时空转动参数a的值，a与r同量纲。
%     返回给定B-L坐标对应的球坐标r、theta值，r_sphere与theta_sphere均为与r_bl、theta_bl同长度的数组。
%     球坐标r = 0处的theta_sphere值始终返回为(0.5 * pi)。
%     只保证对非负实值的r_bl、实值的theta_bl以及实值的a给出正确的返回结果。
    r_sphere = sqrt(r_bl.^2.0 + (a .* sin(theta_bl)).^2.0);
    theta_sphere = acos(r_bl .* cos(theta_bl) ./ r_sphere);
    theta_sphere(isnan(theta_sphere) | isinf(theta_sphere)) = 0.5 * pi;
end
