function [r_extreme, ex_property] = effp_extreme(r_values, is_upper, is_timelike, L, M, a)
% ��Kerrʱ���ⲿ����棨B-L����theta = pi/2������ʱ������˶����ϣ��£���Ч��V+��V-����ʵֵ���ڸ���λ�������еģ����ƣ���ֵ�㡣
% ʹ����Ȼ��λ�ƣ�c = G = 1��
%
% function [r_extreme, ex_property] = effp_extreme(r_values, is_upper, is_timelike, L, M, a)
%     Ҫ����ڼ����Ӧ��Ч��ʵֵ�ε�eff_potential_r������
%     r_valuesָ��λ�����飬ΪB-L����r��ֵ��Ӧ������3��Ԫ�أ�
%     is_upperָ�������Ƿ�Ϊ����Ч�ƣ�����Ϊ����Ч�ƣ�
%     is_timelikeָ�������Ƿ�Ϊ��ʱ��Ч�ƣ�����Ϊ�����Ч�ƣ�
%     Lָ���˶�����B-L����phi���غ�����ֵ��
%     M��a�ֱ�ָ������Դ��������ת������ֵ��a��Mͬ���١�
%     ������ָ������Ч��ʵֵ����ָ��λ�������еĽ��Ƽ�ֵ��r_extreme���Լ���Ӧ��ֵ��ļ�ֵ����ex_property��
%         1Ϊ����ֵ��-1Ϊ��Сֵ��0Ϊ�Ǽ�ֵ�ĵ�������ȷ����
%     ���߾�Ϊ��������ÿ�ж�Ӧһ����ֵ�㡣
    effp_vals = eff_potential_r(r_values(:), is_upper, is_timelike, L, M, a, false); % ӦΪһ������
    val_comps = [sign(effp_vals(2: 1: (end - 1)) - effp_vals(1: 1: (end - 2))), ...
                 sign(effp_vals(2: 1: (end - 1)) - effp_vals(3: 1: end))];
    ex_comp_indexes = find(val_comps(:, 1) == val_comps(:, 2));
    % ע��r_values��effp_vals��Ӧ��val_compsʱ����β����ȥ��һ��Ԫ�ص�
    r_extreme = r_values(ex_comp_indexes + 1);
    ex_property = val_comps(ex_comp_indexes, 1);
end
