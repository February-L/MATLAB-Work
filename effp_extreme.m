function [r_extreme, ex_property] = effp_extreme(r_values, is_upper, is_timelike, L, M, a)
% 求Kerr时空外部赤道面（B-L坐标theta = pi/2）上类时或类光运动的上（下）有效势V+（V-）的实值段在给定位置数组中的（近似）极值点。
% 使用自然单位制：c = G = 1。
%
% function [r_extreme, ex_property] = effp_extreme(r_values, is_upper, is_timelike, L, M, a)
%     要求存在计算对应有效势实值段的eff_potential_r函数。
%     r_values指定位置数组，为B-L坐标r的值，应至少有3个元素；
%     is_upper指定所求是否为上有效势，否则为下有效势；
%     is_timelike指定所求是否为类时有效势，否则为类光有效势；
%     L指定运动关于B-L坐标phi的守恒量的值；
%     M、a分别指定引力源的质量与转动参数值，a与M同量纲。
%     返回所指定的有效势实值段在指定位置数组中的近似极值点r_extreme，以及对应极值点的极值属性ex_property：
%         1为极大值，-1为极小值，0为非极值的导数零点或不确定。
%     两者均为列向量，每行对应一个极值点。
    effp_vals = eff_potential_r(r_values(:), is_upper, is_timelike, L, M, a, false); % 应为一列向量
    val_comps = [sign(effp_vals(2: 1: (end - 1)) - effp_vals(1: 1: (end - 2))), ...
                 sign(effp_vals(2: 1: (end - 1)) - effp_vals(3: 1: end))];
    ex_comp_indexes = find(val_comps(:, 1) == val_comps(:, 2));
    % 注意r_values或effp_vals对应至val_comps时是首尾各截去了一个元素的
    r_extreme = r_values(ex_comp_indexes + 1);
    ex_property = val_comps(ex_comp_indexes, 1);
end
