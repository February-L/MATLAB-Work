function [r_extreme, ex_property, adopted_step] = effp_extreme_s(r_limits, step_values, is_upper, is_timelike, L, M, a)
% 按给定步长序列扫描Kerr时空外部赤道面（B-L坐标theta = pi/2）上类时或类光运动的上（下）有效势V+（V-）的实值段在给定区段内的（近似）极值点。
% 后缀s代表scan。
% 使用自然单位制：c = G = 1。
%
% function [r_extreme, ex_property, adopted_step] = effp_extreme_s(r_limits, step_values, is_upper, is_timelike, L, M, a)
%     要求存在计算对应有效势实值段的eff_potential_r函数。
%     r_limits指定求极值点的区段，应为一2元素数组，第1、2个元素分别指定区段的起点与终点；
%     step_values以数组的形式指定步长序列，若以当前步长扫描出极值点，则缩短扫描区段于极值点附近以下一步长继续扫描，
%     若以当前步长求得某一区段内存在附近有效势值变化不高于double类型精度的极值点，则对该区段仍取上一步长的结果（若为初始区段则舍弃这样的极值点），
%     如此直至步长序列结束，而步长序列应为同号且绝对值递减的，否则扫描过程的正确性不能保证；
%     is_upper指定所求是否为上有效势，否则为下有效势；
%     is_timelike指定所求是否为类时有效势，否则为类光有效势；
%     L指定运动关于B-L坐标phi的守恒量的值；
%     M、a分别指定引力源的质量与转动参数值，a与M同量纲。
%     返回所指定的有效势实值段在指定区段中的近似极值点r_extreme，对应极值点的极值属性ex_property：
%         1为极大值，-1为极小值，
%     以及近似极值点对应的最终采用的扫描步长adopted_step。三者均为列向量，每行对应一个极值点。
    eps_factor = 2.0; % 决定考虑为计算误差的范围的长度的因子，只有整数与半整数意义
    step_count = length(step_values);
    ex_count = 1; % 当前扫描出的极值点个数，也用于r_lower及r_upper的元素个数记忆，故此处先设为1
    r_lower = r_limits(1); r_upper = r_limits(2);
    r_extreme = []; ex_property = []; adopted_step = [];
    for i_step = 1: 1: step_count
        r_ex_tmp = [];
        ex_prop_tmp = [];
        ad_step_tmp = [];
        for i_ex = 1: 1: ex_count
            r_values = (r_lower(i_ex): step_values(i_step): r_upper(i_ex))';
            if (length(r_values) < 3)
                break;
            end
            effp_vals = eff_potential_r(r_values, is_upper, is_timelike, L, M, a, false); % 应为一列向量
            val_comps = [sign(effp_vals(2: 1: (end - 1)) - effp_vals(1: 1: (end - 2))), ...
                         sign(effp_vals(2: 1: (end - 1)) - effp_vals(3: 1: end))];
            ex_comp_indexes = find((val_comps(:, 1) == val_comps(:, 2)) & (val_comps(:, 1) ~= 0.0));
            % 注意r_values或effp_vals对应至val_comps时是首尾各截去了一个元素的
            if (all( ...
                (abs(effp_vals(ex_comp_indexes + 1) - effp_vals(ex_comp_indexes)) > ...
                    (eps_factor * (eps(effp_vals(ex_comp_indexes)) + eps(effp_vals(ex_comp_indexes + 1))))) ...
                & (abs(effp_vals(ex_comp_indexes + 2) - effp_vals(ex_comp_indexes + 1)) > ...
                    (eps_factor * (eps(effp_vals(ex_comp_indexes + 1)) + eps(effp_vals(ex_comp_indexes + 2))))) ...
            ) && (~isempty(ex_comp_indexes)))
                r_ex_tmp = [r_ex_tmp; r_values(ex_comp_indexes + 1)]; %#ok<AGROW>
                ex_prop_tmp = [ex_prop_tmp; val_comps(ex_comp_indexes, 1)]; %#ok<AGROW>
                ad_step_tmp = [ad_step_tmp; step_values(i_step) + zeros(length(ex_comp_indexes), 1)]; %#ok<AGROW>
            elseif (~isempty(r_extreme))
                r_ex_tmp = [r_ex_tmp; r_extreme(i_ex)]; %#ok<AGROW>
                ex_prop_tmp = [ex_prop_tmp; ex_property(i_ex)]; %#ok<AGROW>
                ad_step_tmp = [ad_step_tmp; adopted_step(i_ex)]; %#ok<AGROW>
            elseif (i_step == 1)
                good_index_indexes = find((abs(effp_vals(ex_comp_indexes + 1) - effp_vals(ex_comp_indexes)) > ...
                                       (eps_factor * (eps(effp_vals(ex_comp_indexes)) + eps(effp_vals(ex_comp_indexes + 1))))) ...
                                   & (abs(effp_vals(ex_comp_indexes + 2) - effp_vals(ex_comp_indexes + 1)) > ...
                                       (eps_factor * (eps(effp_vals(ex_comp_indexes + 1)) + eps(effp_vals(ex_comp_indexes + 2))))));
                r_ex_tmp = r_values(ex_comp_indexes(good_index_indexes) + 1);
                ex_prop_tmp = val_comps(ex_comp_indexes(good_index_indexes), 1);
                ad_step_tmp = step_values(i_step) + zeros(length(good_index_indexes), 1);
            end
        end
        r_extreme = r_ex_tmp;
        ex_property = ex_prop_tmp;
        adopted_step = ad_step_tmp;
        ex_count = length(r_extreme);
        if (ex_count == 0)
            break;
        end
        r_lower = r_extreme - step_values(i_step);
        r_upper = r_extreme + step_values(i_step);
        if ((r_lower(1) - r_limits(1)) * sign(step_values(i_step)) < 0.0)
            r_lower(1) = r_limits(1);
        end
        if ((r_upper(end) - r_limits(2)) * sign(step_values(i_step)) > 0.0)
            r_upper(end) = r_limits(2);
        end
    end
end
