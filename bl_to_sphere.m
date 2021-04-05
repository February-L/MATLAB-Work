function [r_sphere, theta_sphere] = bl_to_sphere(r_bl, theta_bl, a)
% ��Kerrʱ��ת������Ϊa��B-L����ת��Ϊ�����꣨���߲�ͬ�Ľ�Ϊr��theta����Ķ��壩��
% ʹ����Ȼ��λ�ƣ�c = G = 1��
%
% function [r_sphere, theta_sphere] = bl_to_sphere(r_bl, theta_bl, a)
%     r_bl��theta_bl�ֱ�ָ��B-L�����r��thetaֵ��ӦΪ������ͬ���ȵ����飻
%     aָ��B-L������ص�Kerrʱ��ת������a��ֵ��a��rͬ���١�
%     ���ظ���B-L�����Ӧ��������r��thetaֵ��r_sphere��theta_sphere��Ϊ��r_bl��theta_blͬ���ȵ����顣
%     ������r = 0����theta_sphereֵʼ�շ���Ϊ(0.5 * pi)��
%     ֻ��֤�ԷǸ�ʵֵ��r_bl��ʵֵ��theta_bl�Լ�ʵֵ��a������ȷ�ķ��ؽ����
    r_sphere = sqrt(r_bl.^2.0 + (a .* sin(theta_bl)).^2.0);
    theta_sphere = acos(r_bl .* cos(theta_bl) ./ r_sphere);
    theta_sphere(isnan(theta_sphere) | isinf(theta_sphere)) = 0.5 * pi;
end
