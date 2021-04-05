function [r_bl, theta_bl] = sphere_to_bl(r_sphere, theta_sphere, a)
% ��������ת��ΪKerrʱ��ת������Ϊa��B-L���꣨���߲�ͬ�Ľ�Ϊr��theta����Ķ��壩��
% ʹ����Ȼ��λ�ƣ�c = G = 1��
%
% function [r_bl, theta_bl] = sphere_to_bl(r_sphere, theta_sphere, a)
%     r_sphere��theta_sphere�ֱ�ָ���������r��thetaֵ��ӦΪ������ͬ���ȵ����飻
%     aָ��B-L������ص�Kerrʱ��ת������a��ֵ��a��rͬ���١�
%     ���ظ����������Ӧ��B-L����r��thetaֵ��r_bl��theta_bl��Ϊ��r_sphere��theta_sphereͬ���ȵ����顣
%     B-L����r = 0����theta_bl����ֵ��������[0, pi/2]�ڡ�
%     ֻ��֤�ԷǸ�ʵֵ��r_sphere��ʵֵ��theta_sphere�Լ�����ʵֵ��a������ȷ�ķ��ؽ����
    r_bl = sqrt(0.5 .* ( ...
               r_sphere.^2.0 - a.^2.0 + ...
               sqrt((r_sphere.^2.0 - a.^2.0).^2.0 + 4.0 .* (a .* r_sphere .* cos(theta_sphere)).^2.0) ...
           ));
    theta_bl = acos(r_sphere .* cos(theta_sphere) ./ r_bl);
    sp_indexes = find(isnan(theta_bl) | isinf(theta_bl));
    theta_bl(sp_indexes) = asin(r_sphere(sp_indexes) ./ abs(a));
end
