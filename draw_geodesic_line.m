clear;
% ���㲢����Kerrʱ���ⲿ����棨B-L����theta = pi/2���ϵ���ʱ���������
% ���³��ֵ�r������δ�ر�˵������ΪB-L�����е�r����
% ��Ȼ��λ�ƣ�c = G = 1
% 
% ��E��Lͬʱ���ţ���dr/ds���䡢dphi/ds���š�dt/ds���ţ�
% ��E��aͬʱ���ţ���dr/ds���䡢dphi/ds���䡢dt/ds���ţ�
% ��a��Lͬʱ���ţ���dr/ds���䡢dphi/ds���š�dt/ds���䡣 

M = 1.0; % ����Դ����
a = 0.8; % ת������
motion_timelike = false; % �Ƿ�Ϊ��ʱ�˶���������Ϊ����˶�
const_L = 5.0; % �˶��������ת����phi���غ���L��ע����a�ϴ�ʱ���������˶����ž޴������

r_horizons = calc_horizon(M, a); % �ӽ�λ��
r_static_surfs = [0.0; 2.0 * M]; % ��ֹ����λ��
V_upper_func = @(r) (eff_potential(r, true, motion_timelike, const_L, M, a, false)); % ����Ч��V+
V_lower_func = @(r) (eff_potential(r, false, motion_timelike, const_L, M, a, false)); % ����Ч��V-
% �������B-L����r���˶�����Ϊ(dr/ds)^2 = g_phiphi * (E - V+(r)) * (E - V-(r)) / r^2, ���жȹ�g_phiphi = r^2 + a^2 + 2 * M * a^2 / r
[r_ex_a, ex_p_a, ad_step_a] = effp_extreme_s([0.0, 100.0], exp((-5.0):(-1.0):(-20.0)), true, motion_timelike, const_L, M, a);
[r_ex_b, ex_p_b, ad_step_b] = effp_extreme_s([0.0, 100.0], exp((-5.0):(-1.0):(-20.0)), false, motion_timelike, const_L, M, a);
% ��Ч�Ƹ������ƣ���ֵ��λ�á���ֵ���ԡ�λ�þ��ȣ���׺a��b�ֱ�ָ�ϡ�����Ч��

draw_with_spherical_r = false; % ������Ϊ�棬��ʹ��������r����B-L����r��ͼ��ע�����ཫ��Ӧ�ظı��ͼ�õ���r������ر�����ֵ

const_E = 1.6; % �˶��������ʱ��t���غ���E��ע����a�ϴ�ʱ���������˶����ž޴������
init_r = r_ex_a(1); % ��ʼrֵ��Ĭ��Ϊ�Ǹ�ֵ
init_dr_direction = '-'; % ��ʼdr/ds����������'+'Ϊ����'-'�����������ַ�Ϊ��
init_phi = 0.0; % ��ʼphiֵ

max_step_len = 0.01; % ΪODE��������õ����s������ע����˶���Χ�ڵ����������s�ı仯������ѡȡ���ʵ�ֵ
s_end = 300.0; r_axis_end = 0.0; % �ֱ�Ϊ�����㵽�����߲���s������ʱ�˶�Ϊ�߳���Ϊ����ʱֹͣ��ͼ��չʾ��B-L����rΪ���٣�������0ʱ���Ӵ˲�����
global phi_end r_end; % �ֱ�Ϊ����phi��rΪĳֵʱֹͣ���㣬��ΪNaN��ֹͣ���ڽӽ��滷ʱ�������ǻ�ֹͣ
phi_end = 100.0 * pi; r_end = 100.0;

draw_multiple_lines = false; % ������Ϊ�棬����õ�ǰ���ڵ���param_names��һ���г����Ƶ������ǰline_count��Ԫ���滻��Ӧ��������һϵ�в����
param_names = ["timelike_vals", "L_vals", "E_vals", "r0_vals", "direction_chars", "max_step_vals";
               "motion_timelike", "const_L", "const_E", "init_r", "init_dr_direction", "max_step_len"];
% ���ַ�������ĵ�һ��ΪӦΪ�������������ڶ���Ϊ����Ԫ��Ҫ��Ӧ�滻�Ĳ����������ƣ�������Ӧ�ַ���������չ
line_count = 2; % ����һϵ�в����ʱ����ߵľ�������
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
polarplot(circle_phi_vals, zeros(1, circle_phi_count) + r_static_surfs(1), 'r-', 'LineWidth', 1.0); % �滷
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
                   (1.0 + (a / init_r)^2.0 + 2.0 * M * a^2.0 / init_r^3.0));  % ��ʼdr/dsֵ
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
    global C13_2m   C14_few  C15_L; % һϵ���˶������еĳ�����ϣ��ڷ��̺���eq_func�����ڼ�����
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
    global C13_2m   C14_few  C15_L; % һϵ���˶������еĳ������
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
    % ���ӽ���dphi/dsΪ���dphi/drΪ���r^(-1)(r����0)�׵�����޷�����ֵ���������ȥ��
    % ������㹻�����ӽ�ʱ���Ե�phi�ı仯����ȥ���ӽ�����phi������������Խ��ӽ������ν���ͬһͼ��
end

function [gstop, isterminal, direction] = stop_event(~, y)
    global use_phi_end phi_end use_r_end r_end;
    gstop = [y(1) - 1E-3; y(1) - r_end; y(3) - phi_end];
    isterminal = [1; use_r_end; use_phi_end];
    direction = [0; 0; 0];
end
