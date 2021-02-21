function horizon_r = calc_horizon(M, g)
    % function horizon_r = calc_horizon(M, g)
    %     ����Bardeenʱ�յ�����Դ����������M�����򻯲���g�������ӽ�뾶��
    %     ʹ����Ȼ��λ�ƣ�c = G = 1, mu_0 = 4*pi��
    %     �������ӽ�뾶�����ڵ��⣩��ɵ������飨Ԫ����Ϊ�ӽ����������ӽ缴���ؿ����顣
    %     (g / M)^2 < 16/27ʱ�����ӽ磬��g������0ʱΪ˫�ӽ�(g = 0���˻���Schwarzschildʱ�գ�r = 0���γ����)��
    metric_roots = sqrt(roots([1.0; 3.0 * g^2.0 - 4.0 * M^2.0; 3.0 * g^4.0; g^6.0]));
    horizon_r = zeros(1, 2);
    for i = 1: 1: 3
        if ((imag(metric_roots(i)) == 0.0) && (real(metric_roots(i)) > 0.0))
            if (horizon_r(1) == 0.0)
                horizon_r(1) = metric_roots(i);
            else
                horizon_r(2) = metric_roots(i);
            end
        end
    end
    horizon_r = unique(horizon_r(horizon_r ~= 0.0));
    horizon_r = sort(horizon_r, 'ascend');
end
