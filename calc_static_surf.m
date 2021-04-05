function r_static = calc_static_surf(the_theta, M, a)
% ����Kerrʱ�վ�ֹ�����ڸ�����λ����B-L����theta��ʾ����λ�ã���B-L����r��ʾ����
% ʹ����Ȼ��λ�ƣ�c = G = 1��
%
% function r_static = calc_static_surf(the_theta, M, a)
%     the_thetaָ����λ��ΪB-L����theta��ֵ�����飩��
%     M��a�ֱ�ָ������Դ��������ת������ֵ��a��Mͬ���١�
%     ����������ֹ�����ڸ�����λ�����У���B-L����rֵ���ɵ�2�о���r_static�����ھ����ڵ�1�У����÷�λ�ϲ����ھ�ֹ�������ӦֵΪNaN��
    if (M <= 0.0)
        r_static = zeros(2, numel(the_theta)) + NaN;
        return;
    end
    the_theta = (the_theta(:))';
    r_static = [M - sqrt(M^2.0 - (a .* cos(the_theta)).^2.0);
                M + sqrt(M^2.0 - (a .* cos(the_theta)).^2.0)];
    r_static(imag(r_static) ~= 0.0) = NaN;
end
