function [r_extreme, ex_property, adopted_step] = effp_extreme_s(r_limits, step_values, is_upper, is_timelike, L, M, a)
% ��������������ɨ��Kerrʱ���ⲿ����棨B-L����theta = pi/2������ʱ������˶����ϣ��£���Ч��V+��V-����ʵֵ���ڸ��������ڵģ����ƣ���ֵ�㡣
% ��׺s����scan��
% ʹ����Ȼ��λ�ƣ�c = G = 1��
%
% function [r_extreme, ex_property, adopted_step] = effp_extreme_s(r_limits, step_values, is_upper, is_timelike, L, M, a)
%     Ҫ����ڼ����Ӧ��Ч��ʵֵ�ε�eff_potential_r������
%     r_limitsָ����ֵ������Σ�ӦΪһ2Ԫ�����飬��1��2��Ԫ�طֱ�ָ�����ε�������յ㣻
%     step_values���������ʽָ���������У����Ե�ǰ����ɨ�����ֵ�㣬������ɨ�������ڼ�ֵ�㸽������һ��������ɨ�裬
%     ���Ե�ǰ�������ĳһ�����ڴ��ڸ�����Ч��ֵ�仯������double���;��ȵļ�ֵ�㣬��Ը�������ȡ��һ�����Ľ������Ϊ��ʼ���������������ļ�ֵ�㣩��
%     ���ֱ���������н���������������ӦΪͬ���Ҿ���ֵ�ݼ��ģ�����ɨ����̵���ȷ�Բ��ܱ�֤��
%     is_upperָ�������Ƿ�Ϊ����Ч�ƣ�����Ϊ����Ч�ƣ�
%     is_timelikeָ�������Ƿ�Ϊ��ʱ��Ч�ƣ�����Ϊ�����Ч�ƣ�
%     Lָ���˶�����B-L����phi���غ�����ֵ��
%     M��a�ֱ�ָ������Դ��������ת������ֵ��a��Mͬ���١�
%     ������ָ������Ч��ʵֵ����ָ�������еĽ��Ƽ�ֵ��r_extreme����Ӧ��ֵ��ļ�ֵ����ex_property��
%         1Ϊ����ֵ��-1Ϊ��Сֵ��
%     �Լ����Ƽ�ֵ���Ӧ�����ղ��õ�ɨ�貽��adopted_step�����߾�Ϊ��������ÿ�ж�Ӧһ����ֵ�㡣
    eps_factor = 2.0; % ��������Ϊ�������ķ�Χ�ĳ��ȵ����ӣ�ֻ�����������������
    step_count = length(step_values);
    ex_count = 1; % ��ǰɨ����ļ�ֵ�������Ҳ����r_lower��r_upper��Ԫ�ظ������䣬�ʴ˴�����Ϊ1
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
            effp_vals = eff_potential_r(r_values, is_upper, is_timelike, L, M, a, false); % ӦΪһ������
            val_comps = [sign(effp_vals(2: 1: (end - 1)) - effp_vals(1: 1: (end - 2))), ...
                         sign(effp_vals(2: 1: (end - 1)) - effp_vals(3: 1: end))];
            ex_comp_indexes = find((val_comps(:, 1) == val_comps(:, 2)) & (val_comps(:, 1) ~= 0.0));
            % ע��r_values��effp_vals��Ӧ��val_compsʱ����β����ȥ��һ��Ԫ�ص�
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
