clear;
% 计算Bardeen时空类时及类光测地线（轨道面theta = pi/2）
% 自然单位制：c = G = 1, mu_0 = 4*pi
% Bardeen时空非零度规：g00 = -f(r), g11 = 1 / f(r), g22 = r^2, g33 = (r * sin(theta))^2，
% 其中f(r) = 1 - 2 * m * r^2 / (r^2 + g^2)^(3/2)。

global motion_timelike; % 是否为类时运动，否则视为类光运动
motion_timelike = false;
M = 1.0; % 引力源质量，对史瓦西时空而言视界位于r = 2M处
g = 0.7; % 正则化参数，可解释为一非线性磁单极子的磁荷，在(g/M)^2 < 16/27（即g小于约0.7698M）时Bardeen时空存在两个视界
const_L = 3.6; % 物体关于转动角phi的守恒量L，L = r^2 * d(phi)/ds，对类光运动s换为非零仿射参数

r_horizons = calc_horizon(M, g); % （各）视界半径
if (motion_timelike) % 类时运动
    V2_func = @(r) ((1.0 - 2.0 .* M .* r.^2.0 ./ (r.^2.0 + g.^2.0).^1.5) .* (1.0 + (const_L./r).^2.0));
else % 类光运动
    V2_func = @(r) ((1.0 - 2.0 .* M .* r.^2.0 ./ (r.^2.0 + g.^2.0).^1.5) .* (const_L./r).^2.0);
end
% 有效势的平方V^2，(dr/ds)^2 + V^2 = E^2
% 若忽略V^2表达式中左右两因子各自第二项的乘积（在光速c趋于无穷时的更高阶小量）的贡献，则该有效势导出的轨道方程成为牛顿近似下的方程。
r_extremes = calc_V2_extreme(M, g, const_L, motion_timelike); % 有效势（各）极值点位置
%r_vals = 0.1: 0.01: 25.0; V2_vals = V2_func(r_vals);
%figure(1); plot(r_vals, V2_vals, "k-"); xlabel("r / r_S"); ylabel("V^2"); axis([1.0, 25.0, 0.9, ceil(max(V2_vals)*10.0)*0.1]);

const_E2 = 0.55; % 物体关于时间t的守恒量的平方E^2，E = f(r) * dt/ds，对类光运动s换为非零仿射参数
E2_vals = [1.2, 1.6, 2.0, 2.4];
use_E2_vals = false; % 若该项为真，则不采用上面的const_E2而是采用E2_vals的值求一系列解

init_r = r_extremes(1); % phi = 0时的r值
init_dr_direction = -1.0; % phi = 0时dr/d(phi)的正负方向，+1.0为正，-1.0为负
init_dr = init_dr_direction * sqrt(const_E2 - V2_func(init_r)) / const_L * init_r^2.0; % phi = 0时的dr/d(phi)值
global g2 factor_M factor_g factor_L;
g2 = g^2.0; factor_M = 3.0 * M; factor_g =  2.0 * M * g^2.0 / const_L^2.0; factor_L = M / const_L^2.0;
% 若强行取factor_M = 0.0，则方程变为牛顿近似下的轨道方程。

phi_end = 80.0 * pi; r_axis_end = 0.0; % 分别为：计算到phi为多少时停止、图像展示到r为多少（不大于0时无视此参数）
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
    if (motion_timelike) % 类时运动
        du_and_dp = [y(2);
                     2.0 .* y(2).^2.0 ./ y(1) + y(1) + y(1).^5.0 ./ (y(1).^2.0 + g2).^2.5 .* (factor_g - factor_M - factor_L .* y(1).^2.0)];
    else % 类光运动
        du_and_dp = [y(2);
                     2.0 .* y(2).^2.0 ./ y(1) + y(1) - factor_M .* y(1).^5.0 ./ (y(1).^2.0 + g2).^2.5];
    end
end

function [gstop, isterminal, direction] = stop_event(~, y)
    gstop = [y(1) - 1E-8; y(1) - 100.0];
    isterminal = [1; 1];
    direction = [0; 0];
end
