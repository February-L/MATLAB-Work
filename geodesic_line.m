clear;
% ����Bardeenʱ����ʱ��������ߣ������theta = pi/2��
% ��Ȼ��λ�ƣ�c = G = 1, mu_0 = 4*pi
% Bardeenʱ�շ���ȹ棺g00 = -f(r), g11 = 1 / f(r), g22 = r^2, g33 = (r * sin(theta))^2��
% ����f(r) = 1 - 2 * m * r^2 / (r^2 + g^2)^(3/2)��

global motion_timelike; % �Ƿ�Ϊ��ʱ�˶���������Ϊ����˶�
motion_timelike = false;
M = 1.0; % ����Դ��������ʷ����ʱ�ն����ӽ�λ��r = 2M��
g = 0.7; % ���򻯲������ɽ���Ϊһ�����Դŵ����ӵĴźɣ���(g/M)^2 < 16/27����gС��Լ0.7698M��ʱBardeenʱ�մ��������ӽ�
const_L = 3.6; % �������ת����phi���غ���L��L = r^2 * d(phi)/ds��������˶�s��Ϊ����������

r_horizons = calc_horizon(M, g); % �������ӽ�뾶
if (motion_timelike) % ��ʱ�˶�
    V2_func = @(r) ((1.0 - 2.0 .* M .* r.^2.0 ./ (r.^2.0 + g.^2.0).^1.5) .* (1.0 + (const_L./r).^2.0));
else % ����˶�
    V2_func = @(r) ((1.0 - 2.0 .* M .* r.^2.0 ./ (r.^2.0 + g.^2.0).^1.5) .* (const_L./r).^2.0);
end
% ��Ч�Ƶ�ƽ��V^2��(dr/ds)^2 + V^2 = E^2
% ������V^2���ʽ�����������Ӹ��Եڶ���ĳ˻����ڹ���c��������ʱ�ĸ��߽�С�����Ĺ��ף������Ч�Ƶ����Ĺ�����̳�Ϊţ�ٽ����µķ��̡�
r_extremes = calc_V2_extreme(M, g, const_L, motion_timelike); % ��Ч�ƣ�������ֵ��λ��
%r_vals = 0.1: 0.01: 25.0; V2_vals = V2_func(r_vals);
%figure(1); plot(r_vals, V2_vals, "k-"); xlabel("r / r_S"); ylabel("V^2"); axis([1.0, 25.0, 0.9, ceil(max(V2_vals)*10.0)*0.1]);

const_E2 = 0.55; % �������ʱ��t���غ�����ƽ��E^2��E = f(r) * dt/ds��������˶�s��Ϊ����������
E2_vals = [1.2, 1.6, 2.0, 2.4];
use_E2_vals = false; % ������Ϊ�棬�򲻲��������const_E2���ǲ���E2_vals��ֵ��һϵ�н�

init_r = r_extremes(1); % phi = 0ʱ��rֵ
init_dr_direction = -1.0; % phi = 0ʱdr/d(phi)����������+1.0Ϊ����-1.0Ϊ��
init_dr = init_dr_direction * sqrt(const_E2 - V2_func(init_r)) / const_L * init_r^2.0; % phi = 0ʱ��dr/d(phi)ֵ
global g2 factor_M factor_g factor_L;
g2 = g^2.0; factor_M = 3.0 * M; factor_g =  2.0 * M * g^2.0 / const_L^2.0; factor_L = M / const_L^2.0;
% ��ǿ��ȡfactor_M = 0.0���򷽳̱�Ϊţ�ٽ����µĹ�����̡�

phi_end = 80.0 * pi; r_axis_end = 0.0; % �ֱ�Ϊ�����㵽phiΪ����ʱֹͣ��ͼ��չʾ��rΪ���٣�������0ʱ���Ӵ˲�����
max_step_len = 0.01 * pi;
the_canvas = figure(2);
circle_phi_vals = 0.0: max_step_len: (2.0*pi); circle_phi_count = length(circle_phi_vals);
polarplot(circle_phi_vals, zeros(1, circle_phi_count) + r_horizons(1), "r:");
hold on;
polarplot(circle_phi_vals, zeros(1, circle_phi_count) + r_horizons(2), "b:");

ode_options = odeset('MaxStep', max_step_len, 'AbsTol', 1E-8, 'events', @stop_event);
if (~use_E2_vals)
    the_solution = ode45(@eq_func, [0.0, phi_end], [init_r; init_dr], ode_options);
    polarplot(the_solution.x, the_solution.y(1, :), "k-");
else
    for i = 1: 1: length(E2_vals)
        init_dr = init_dr_direction * sqrt(E2_vals(i) - V2_func(init_r)) / const_L * init_r^2.0;
        the_solution = ode45(@eq_func, [0.0, phi_end], [init_r; init_dr], ode_options);
        polarplot(the_solution.x, the_solution.y(1, :), "k-");
    end
end
if (r_axis_end > 0.0)
    axis_now = axis();
    axis_now(4) = r_axis_end;
    axis(axis_now);
end
hold off;
set(get(the_canvas, 'Children'), 'FontName', 'Times New Roman');

function du_and_dp = eq_func(~, y)
    global motion_timelike g2 factor_M factor_g factor_L;
    % y(1) is r, y(2) is p = dr/d(phi).
    if (motion_timelike) % ��ʱ�˶�
        du_and_dp = [y(2);
                     2.0 .* y(2).^2.0 ./ y(1) + y(1) + y(1).^5.0 ./ (y(1).^2.0 + g2).^2.5 .* (factor_g - factor_M - factor_L .* y(1).^2.0)];
    else % ����˶�
        du_and_dp = [y(2);
                     2.0 .* y(2).^2.0 ./ y(1) + y(1) - factor_M .* y(1).^5.0 ./ (y(1).^2.0 + g2).^2.5];
    end
end

function [gstop, isterminal, direction] = stop_event(~, y)
    gstop = [y(1) - 1E-8; y(1) - 100.0];
    isterminal = [1; 1];
    direction = [0; 0];
end
