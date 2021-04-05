function effp_val = eff_potential_r(r, is_upper, is_timelike, L, M, a, clean_imag_err)
% ����Kerrʱ���ⲿ����棨B-L����theta = pi/2������ʱ������˶����ϣ��£���Ч��V+��V-���ڸ���λ�õ�ֵ��
% ��׺r����real����ֵ����ΪNaN��
% ʹ����Ȼ��λ�ƣ�c = G = 1��
%
% function effp_val = eff_potential_r(r, is_upper, is_timelike, L, M, a, clean_imag_err)
%     rָ��λ�ã�ΪB-L����r��ֵ��
%     is_upperָ�������Ƿ�Ϊ����Ч�ƣ�����Ϊ����Ч�ƣ�
%     is_timelikeָ�������Ƿ�Ϊ��ʱ��Ч�ƣ�����Ϊ�����Ч�ƣ�
%     Lָ���˶�����B-L����phi���غ�����ֵ��
%     M��a�ֱ�ָ������Դ��������ת������ֵ��a��Mͬ���٣�
%     clean_imag_errָ���Ƿ񽫣��ӽ總�����������µ��鲿Ĩȥ��
%     ������ָ������Ч����ָ��λ�ã����飩��ֵ�����飩effp_val����ֵ����ΪNaN��
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
