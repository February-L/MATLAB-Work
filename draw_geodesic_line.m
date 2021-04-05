clear;
% 计算并绘制Kerr时空外部赤道面（B-L坐标theta = pi/2）上的类时及类光测地线
% 以下出现的r坐标若未特别说明，均为B-L坐标中的r坐标
% 自然单位制：c = G = 1
% 
% 若E、L同时反号，则dr/ds不变、dphi/ds反号、dt/ds反号；
% 若E、a同时反号，则dr/ds不变、dphi/ds不变、dt/ds反号；
% 若a、L同时反号，则dr/ds不变、dphi/ds反号、dt/ds不变。 

M = 1.0; % 引力源质量
a = 0.8; % 转动参数
motion_timelike = false; % 是否为类时运动，否则视为类光运动
const_L = 5.0; % 运动物体关于转动角phi的守恒量L，注意在a较大时其正负对运动有着巨大的区别

r_horizons = calc_horizon(M, a); % 视界位置
r_static_surfs = [0.0; 2.0 * M]; % 静止界限位置
V_upper_func = @(r) (eff_potential(r, true, motion_timelike, const_L, M, a, false)); % 上有效势V+
V_lower_func = @(r) (eff_potential(r, false, motion_timelike, const_L, M, a, false)); % 下有效势V-
% 物体关于B-L坐标r的运动方程为(dr/ds)^2 = g_phiphi * (E - V+(r)) * (E - V-(r)) / r^2, 其中度规g_phiphi = r^2 + a^2 + 2 * M * a^2 / r
[r_ex_a, ex_p_a, ad_step_a] = effp_extreme_s([0.0, 100.0], exp((-5.0):(-1.0):(-20.0)), true, motion_timelike, const_L, M, a);
[r_ex_b, ex_p_b, ad_step_b] = effp_extreme_s([0.0, 100.0], exp((-5.0):(-1.0):(-20.0)), false, motion_timelike, const_L, M, a);
% 有效势各（近似）极值点位置、极值属性、位置精度，后缀a、b分别指上、下有效势

draw_with_spherical_r = false; % 若该项为真，则使用球坐标r而非B-L坐标r绘图，注意这亦将对应地改变绘图用到的r坐标相关变量的值

const_E = 1.6; % 运动物体关于时间t的守恒量E，注意在a较大时其正负对运动有着巨大的区别
init_r = r_ex_a(1); % 初始r值，默认为非负值
init_dr_direction = '-'; % 初始dr/ds的正负方向，'+'为正，'-'或其他单个字符为负
init_phi = 0.0; % 初始phi值

max_step_len = 0.01; % 为ODE求解器设置的最大s步长，注意据运动范围内的坐标相对于s的变化快慢来选取合适的值
s_end = 300.0; r_axis_end = 0.0; % 分别为：计算到世界线参数s（对类时运动为线长）为多少时停止、图像展示到B-L坐标r为多少（不大于0时无视此参数）
global phi_end r_end; % 分别为：在phi、r为某值时停止计算，若为NaN则不停止，在接近奇环时计算总是会停止
phi_end = 100.0 * pi; r_end = 100.0;

draw_multiple_lines = false; % 若该项为真，则采用当前存在的在param_names第一行列出名称的数组的前line_count个元素替换对应参数绘制一系列测地线
param_names = ["timelike_vals", "L_vals", "E_vals", "r0_vals", "direction_chars", "max_step_vals";
               "motion_timelike", "const_L", "const_E", "init_r", "init_dr_direction", "max_step_len"];
% 该字符串矩阵的第一行为应为数组命的名，第二行为数组元素要对应替换的参数变量名称，补充相应字符串即可拓展
line_count = 2; % 绘制一系列测地线时测地线的具体条数
%E_vals = [V_upper_func(r_ex_a(3)), V_lower_func(r_ex_b(3))];
%r0_vals = [r_ex_a(3), r_ex_b(3)];

the_canvas = figure(2);
bl_to_sph_func = @(r) (sqrt(r.^2.0 + a^2.0));
circle_phi_vals = 0.0: (0.01 * pi): (2.0*pi); circle_phi_count = length(circle_phi_vals);
if (draw_with_spherical_r)
    r_horizons = bl_to_sph_func(r_horizons);
    r_static_surfs = bl_to_sph_func(r_static_surfs);
end
polarplot(NaN, NaN); hold on;
polarplot(circle_phi_vals, zeros(1, circle_phi_count) + r_horizons(1), 'r--');
polarplot(circle_phi_vals, zeros(1, circle_phi_count) + r_horizons(2), 'b--');
polarplot(circle_phi_vals, zeros(1, circle_phi_count) + r_static_surfs(1), 'r-', 'LineWidth', 1.0); % 奇环
polarplot(circle_phi_vals, zeros(1, circle_phi_count) + r_static_surfs(2), 'g:', 'LineWidth', 1.0);
clear circle_phi_vals circle_phi_count;

global use_phi_end use_r_end;
direction_factor_func = @(drctn_char) (2.0 * double(drctn_char == '+') - 1.0);
if (~draw_multiple_lines)
    loop_count = 1;
else
    loop_count = line_count;
    param_count = size(param_names, 2);
end
for i = 1: 1: loop_count
    if (draw_multiple_lines)
        for i_p = 1: 1: param_count
            if (exist(param_names(1, i_p), 'var'))
                eval(strcat(param_names(2, i_p), "=", param_names(1, i_p), "(i);"));
            end
        end
    end
    init_dr = direction_factor_func(init_dr_direction) * ...
              sqrt((const_E - eff_potential(init_r, true, motion_timelike, const_L, M, a, false)) * ...
                   (const_E - eff_potential(init_r, false, motion_timelike, const_L, M, a, false)) * ...
                   (1.0 + (a / init_r)^2.0 + 2.0 * M * a^2.0 / init_r^3.0));  % 初始dr/ds值
    set_eq_constants(M, a, motion_timelike, const_L, const_E);
    use_phi_end = double(~isnan(phi_end)); use_r_end = double(~isnan(r_end));
    ode_options = odeset('MaxStep', max_step_len, 'InitialStep', max(max_step_len, 1E-3), 'AbsTol', 1E-8, 'events', @stop_event);
    the_solution = ode23s(@eq_func, [0.0, s_end], [init_r; init_dr; init_phi], ode_options);
    if (draw_with_spherical_r)
        the_solution.y(1, :) = bl_to_sph_func(the_solution.y(1, :));
    end
    polarplot(the_solution.y(3, :), the_solution.y(1, :), 'k-');
    fprintf("Geodesic line #%d: M = %G, a = %G, %s, L = %G, E = %.6G.\n", int32(i), M, a, motion_type_str(motion_timelike), const_L, const_E);
end
axis_now = axis();
if (r_axis_end > 0.0)
    if (draw_with_spherical_r)
        r_axis_end = bl_to_sph_func(r_axis_end);
    end
    axis_now(4) = r_axis_end;
    axis(axis_now);
end
hold off;
set(get(the_canvas, 'Children'), 'FontName', 'Times New Roman');

function the_str = motion_type_str(motion_timelike)
    if (motion_timelike)
        the_str = "time-like";
    else
        the_str = "light-like";
    end
end

function set_eq_constants(M, a, motion_timelike, const_L, const_E)
    global C1_a2    C2_3m    C3_E2    C4_2ma2;
    global C5_2mid  C6_L2    C7_2km   C8_k;
    global C9_3mid  C10_many C11_km   C12_4ma2;
    global C13_2m   C14_few  C15_L; % 一系列运动方程中的常数组合，在方程函数eq_func中用于简化内容
    C1_a2 = a^2.0; C2_3m = 3.0*M; C3_E2 = const_E^2.0; C4_2ma2 = 2.0*M*a^2.0;
    C5_2mid = 2.0*M*const_L * (2.0*a*const_E - const_L); C6_L2 = const_L^2.0;
    C8_k = -1.0*double(logical(motion_timelike)); C7_2km = 2.0*C8_k*M;
    C9_3mid = 3.0*M*const_L * (2.0*a*const_E - const_L); C10_many = 2.0*M*a^2.0 * const_L*(a*const_E - const_L);
    C11_km = C8_k*M; C12_4ma2 = 4.0*M*a^2.0; C13_2m = 2.0*M;
    C14_few = 2.0*M * (a*const_E - const_L); C15_L = const_L;
end

function dr_dp_dphi = eq_func(~, y)
    global C1_a2    C2_3m    C3_E2    C4_2ma2;
    global C5_2mid  C6_L2    C7_2km   C8_k;
    global C9_3mid  C10_many C11_km   C12_4ma2;
    global C13_2m   C14_few  C15_L; % 一系列运动方程中的常数组合
    % x is s, y(1) is r, y(2) is p = dr/ds, y(3) is phi
    y1_2 = y(1).^2.0; y1_3 = y(1).^3.0;
    horizon_omega = y1_2 - C13_2m .* y(1) + C1_a2;
    dr_dp_dphi = [y(2);
                 ( ...
                     (C6_L2 .* y1_3 + C9_3mid .* y1_2 + C10_many + C11_km .* ((y1_2 + C1_a2).^2.0 - C12_4ma2 .* y(1))) ...
                         ./ (y1_3 .* (y1_3 + C1_a2 .* y(1) + C4_2ma2)) ...
                     - (C3_E2 - (C5_2mid + C6_L2 .* y(1) + C7_2km .* (y1_2 + C1_a2)) ./ (y1_3 + C1_a2 .* y(1) + C4_2ma2) + C8_k) ...
                         .* C1_a2 ./ y1_3 .* (1.0 + C2_3m ./ y(1)) ...
                 );
                 (C14_few ./ y(1) + C15_L) ./ horizon_omega];
    if (abs(horizon_omega) <= 1E-7)
        dr_dp_dphi(3) = 0.0;
    end
    % 在视界上dphi/ds为无穷，dphi/dr为最低r^(-1)(r趋于0)阶的无穷，无法将数值计算进行下去，
    % 因此在足够靠近视界时忽略掉phi的变化（即去掉视界内外phi间相差的无穷），以将视界两侧衔接于同一图上
end

function [gstop, isterminal, direction] = stop_event(~, y)
    global use_phi_end phi_end use_r_end r_end;
    gstop = [y(1) - 1E-3; y(1) - r_end; y(3) - phi_end];
    isterminal = [1; use_r_end; use_phi_end];
    direction = [0; 0; 0];
end
