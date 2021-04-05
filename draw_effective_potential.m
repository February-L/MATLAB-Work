clear;
format compact; format long g;
% ����Kerrʱ���ⲿ����棨B-L����theta = pi/2����Ч��
% ��Ȼ��λ�ƣ�c = G = 1
% V+(r, L) = -V-(r, -L)��V-(r, L) = -V+(r, -L)����L��Ϊaͬ��

motion_timelike = true; % �Ƿ�Ϊ��ʱ�˶���������Ϊ����˶�
M = 1.0; % ����Դ����
a = 0.8; % ת������
const_L = 5.0; % �˶��������ת����phi���غ���L
var_choice = 'L'; % 'L' - �仯L�̶�M��a��'a' - �仯a�̶�M��L����δ���룩'M' - �仯M�̶�a��L
var_vals = []; % var_choice��var_vals��Ӧ����δ���ƣ��ݲ�ʹ��
var_count = length(var_vals); % ��var_countΪ0�򲻶�M��a��L��var_vals���¸�ֵ
drawing_parts = logical([1, 1]); % �����Ƿ�����ϡ�����Ч�ƣ���һ������Ԫ�طֱ��Ӧ�ϡ�����Ч��
colour_diff = true; % �����Ƿ�Ϊ�ϡ�����Ч�����ֻ�����ɫ
mark_showed_extreme = false; % �����Ƿ��ʾ��������ͼ�еļ�ֵ��
also_farther_extreme = false; % �ڱ�ʾ��ֵ��ʱ�����Ƿ�Ҳ��δ������ͼ�е��ڸ�Զ���ļ�ֵ����
mark_horizon = false; % �����Ƿ��ʾ���ӽ�

r_lower = 0.0; r_upper = 10.0;
r_step_len = 0.01;
V_lower = []; V_upper = []; % V�������������ƣ���ֵΪ���򲻲�ȡ��Ӧ����

horizon_r = calc_horizon(M, a) %#ok<NOPTS>
h1_tick = horizon_r(1); h2_tick = horizon_r(2);
r_init_ticks = r_lower: r_step_len: r_upper;
inner_tick_i = find(r_init_ticks <= h1_tick);
outer_tick_i = find(r_init_ticks >= h2_tick);
mid_tick_i = find(~((r_init_ticks <= h1_tick) | (r_init_ticks >= h2_tick)));
h1_tick(isempty(inner_tick_i) || isempty(mid_tick_i)) = []; h2_tick(isempty(outer_tick_i) || isempty(mid_tick_i)) = [];
r_ticks = [r_init_ticks(inner_tick_i), h1_tick, r_init_ticks(mid_tick_i), h2_tick, r_init_ticks(outer_tick_i)];
clear r_init_ticks h1_tick h2_tick inner_tick_i outer_tick_i mid_tick_i;

[r_ex_a, ex_p_a, ad_step_a] = effp_extreme_s([0.0, 100.0], exp((-5.0):(-1.0):(-20.0)), true, motion_timelike, const_L, M, a);
[r_ex_b, ex_p_b, ad_step_b] = effp_extreme_s([0.0, 100.0], exp((-5.0):(-1.0):(-20.0)), false, motion_timelike, const_L, M, a);
upper_extreme = [r_ex_a, ex_p_a, ad_step_a] %#ok<NASGU,NOPTS>
lower_extreme = [r_ex_b, ex_p_b, ad_step_b] %#ok<NASGU,NOPTS>
upper_extreme_values = eff_potential(r_ex_a, true, motion_timelike, const_L, M, a, false)' %#ok<NOPTS>
lower_extreme_values = eff_potential(r_ex_b, false, motion_timelike, const_L, M, a, false)' %#ok<NOPTS>
clear upper_extreme lower_extreme;

upper_rgb = [0.0, 0.45, 0.74] .* double(colour_diff);
lower_rgb = [0.85, 0.33, 0.1] .* double(colour_diff);
var_names = cell(1, var_count);
loop_time = var_count + (var_count == 0);
effective_choice = lower(var_choice);
the_canvas = figure();
hold on;
for i = 1: 1: loop_time
    if (var_count ~= 0)
        switch (effective_choice)
            case 'l'
                const_L = var_vals(i);
            case 'a'
                a = var_vals(i);
            otherwise
                fprintf("Error: unavailable variable choice.\n");
                return;
        end
        var_names{i} = sprintf("\\fontsize{10}\\it%c\\rm = %g", var_choice, var_vals(i));
        plot(r_ticks, eff_potential_r(r_ticks, true, motion_timelike, const_L, M, a, true), ...
             r_ticks, eff_potential_r(r_ticks, false, motion_timelike, const_L, M, a, true));
    else
        if (drawing_parts(1))
            plot(r_ticks, eff_potential_r(r_ticks, true, motion_timelike, const_L, M, a, true), '-', 'Color', upper_rgb);
        end
        if (drawing_parts(2))
            plot(r_ticks, eff_potential_r(r_ticks, false, motion_timelike, const_L, M, a, true), '-', 'Color', lower_rgb);
        end
    end
end
axis_now = axis();
axis_now([1, 2]) = [r_lower, r_upper];
if (~isempty(V_lower))
    axis_now(3) = V_lower;
end
if (~isempty(V_upper))
    axis_now(4) = V_upper;
end
axis(axis_now);
if (mark_showed_extreme)
    r_ex = [r_ex_b; r_ex_a]';
    V_ex = [lower_extreme_values, upper_extreme_values];
    of_upper = [logical(0.0 .* r_ex_b); ~logical(0.0 .* r_ex_a)]';
    ex_count = length(r_ex);
    for i = 1: 1: ex_count
        if ((r_ex(i) > axis_now(1)) && ((r_ex(i) < axis_now(2)) || also_farther_extreme))
            if (of_upper(i))
                line_rgb = upper_rgb;
            else
                line_rgb = lower_rgb;
            end
            if ((of_upper(i) && drawing_parts(1)) || ((~of_upper(i)) && drawing_parts(2)))
                plot([axis_now(1), r_ex(i)], [V_ex(i), V_ex(i)], '--', 'Color', line_rgb, 'LineWidth', 0.25);
            end
        end
    end
end
if (mark_horizon)
    for i = 1: 1: 2
        if ((horizon_r(i) > axis_now(1)) && (horizon_r(i) < axis_now(2)))
            plot([horizon_r(i), horizon_r(i)], [axis_now(3), eff_potential(horizon_r(i), true, motion_timelike, const_L, M, a, true)], ...
                 'k:', 'LineWidth', 0.75);
        end
    end
end
hold off;
if (var_count ~= 0)
    legend(var_names);
end
xlabel("\it{r}"); ylabel("\it{V}\rm_{\pm}", 'Rotation', 0.0);
set(get(the_canvas, 'Children'), 'FontName', 'Times New Roman');
set(get(the_canvas, 'Children'), 'Box', true);
