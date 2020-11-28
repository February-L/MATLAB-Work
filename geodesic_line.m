clear;
% ����ʷ����ʱ����ʱ����ߣ������theta = pi/2��

r_S = 1.0; % ʷ�����뾶r_S = 2GM/c^2
const_L = 2.3; % �������ת����phi���غ���L��L = r^2 * d(phi)/ds��L^2����3(r_S^2)ʱ�������ȶ�Բ���r = 3 * r_S

V2_func = @(r) ((1.0 - r_S./r) .* (1.0 + (const_L./r).^2.0)); % ��Ч�Ƶ�ƽ��V^2��(dr/ds)^2 + V^2 = E^2
% ������V^2���ʽ�е�[(r_S * L^2) / r^3]��ڹ���c��������ʱ�ĸ��߽�С�����������Ч�Ƶ����Ĺ�����̳�Ϊţ�������µ�Binet���̡�
lower_r0 = const_L.^2.0 / r_S - const_L * sqrt((const_L/r_S).^2.0 - 3.0); % ��Ч�Ƽ���ֵ�㣬����
upper_r0 = const_L.^2.0 / r_S + const_L * sqrt((const_L/r_S).^2.0 - 3.0); % ��Ч�Ƽ�Сֵ�㣬����
%r_vals = 1.0: 0.1: 25.0; V2_vals = V2_func(r_vals);
%figure(1); plot(r_vals, V2_vals, "k-"); xlabel("r / r_S"); ylabel("V^2"); axis([1.0, 25.0, 0.9, ceil(max(V2_vals)*10.0)*0.1]);

const_E2 = V2_func(5); % �������ʱ��ct���غ���E^2��E = (1 - r_S/r) * d(ct)/ds��E >= 1�ſ��ܴ��ӽ����Զ������������Զ
E2_vals = [1.171, 1.175, 1.185, 1.25];
use_E2_vals = false; % ������Ϊ�棬�򲻲��������const_E2���ǲ���E2_vals��ֵ��һϵ�н�

init_r = upper_r0; % phi = 0ʱ��rֵ
init_dr_direction = -1.0; % phi = 0ʱdr/d(phi)����������+1.0Ϊ����-1.0Ϊ��
init_dr = init_dr_direction * sqrt(const_E2 - V2_func(init_r)) / const_L * init_r^2.0; % phi = 0ʱ��dr/d(phi)ֵ
global factor_r factor_L;
factor_r = 1.5 * r_S; factor_L = 0.5 * r_S/(const_L*const_L);
% ��ǿ��ȡfactor_r = 0.0���򷽳̱�Ϊţ�������µó���Binet���̡�

phi_end = 100.0 * pi; r_axis_end = 0.0; % �ֱ�Ϊ�����㵽phiΪ����ʱֹͣ��ͼ��չʾ��rΪ���٣�С��r_S��0ʱ���Ӵ˲�����
max_step_len = 0.01 * pi;
figure(2);
circle_phi_vals = 0.0: max_step_len: (2.0*pi); circle_phi_count = length(circle_phi_vals);
polarplot(circle_phi_vals, zeros(1, circle_phi_count) + r_S, "r-");
hold on;
polarplot(circle_phi_vals, zeros(1, circle_phi_count) + lower_r0, "b:");

ode_options = odeset('MaxStep', max_step_len, 'AbsTol', 1e-8, 'events', @stop_event);
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
if ((r_axis_end >= r_S) && (r_axis_end > 0.0))
    axis_now = axis();
    axis_now(4) = r_axis_end;
    axis(axis_now);
end
hold off;

function du_and_dp = eq_func(~, y)
    global factor_r factor_L;
    % y(1) is r, y(2) is p = dr/d(phi).
    du_and_dp = [y(2);
                 -factor_r + y(1) - factor_L .* y(1).^2.0 + 2.0 .* y(2).^2.0 ./ y(1)];
end

function [gstop, isterminal, direction] = stop_event(~, y)
    gstop = [y(1) - 1.0; y(1) - 100.0];
    isterminal = [1; 1];
    direction = [0; 0];
end
