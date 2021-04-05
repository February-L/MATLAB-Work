function r_horizon = calc_horizon(M, a)
% ����Kerrʱ���ӽ��λ�ã���B-L����r��ʾ����
% ʹ����Ȼ��λ�ƣ�c = G = 1��
%
% function r_horizon = calc_horizon(M, a)
%     M��a�ֱ�ָ������Դ��������ת������ֵ��a��Mͬ���١�
%     ���������ӽ��B-L����rֵ���ɵ�2Ԫ������r_horizon�������ӽ���ǰ�����ӽ粻�������ӦֵΪNaN��
%     ��|a| = M��ʹ���ӽ���r = M����ȫ�غ�ʱ�����ص�r_horizon��Ϊ[M; NaN]��
    a = abs(a);
    if (a < M)
        r_horizon = M + sqrt(M^2.0 - a^2.0) .* [-1.0; 1.0];
    elseif ((a == M) && (M ~= 0.0))
        r_horizon = [M; NaN];
    else
        r_horizon = [NaN; NaN];
    end
end
