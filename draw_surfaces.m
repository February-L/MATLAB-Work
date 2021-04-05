clear;
% ��B-L��������������Kerrʱ�ճ�������ת��������ϵ��ӽ硢��ֹ���޼��滷����
% ��Ȼ��λ�ƣ�c = G = 1

M = 1.0; % ����Դ����
a = 1.2; % ת������
drawing_surfaces = logical([1, 1, 1]); % �����Ƿ���ƶ�Ӧ���棬��һ����������Ԫ�طֱ��Ӧ�ӽ硢��ֹ���ޡ��滷
is_equatorial = false; % ��ȡ���ƽ����Ƿ�Ϊ����棬����Ϊ��ת����Ľ���
is_bl = false; % ��ȡ���������Ƿ�ΪB-L���꣬����Ϊ������
solid_stat_surf = true; % �Ƿ񽫾�ֹ���޻�Ϊʵ�ߣ�����Ϊ�����

r_axis_end = 2.5; % ָ��ͼ��չʾ��rΪ���٣�������0ʱ���Ӵ˲�����
the_canvas = figure();
polarplot(NaN, NaN); hold on;
angle_step = 0.005 * pi;
init_angle_series = 0.0: angle_step: pi;
angle_count = size(init_angle_series, 2);
if (drawing_surfaces(1))
    angle_vals = ones(2, 1) * init_angle_series;
    r_vals = calc_horizon(M, a) * ones(1, angle_count);
    if (~is_bl)
        if (is_equatorial)
            r_vals = bl_to_sphere(r_vals, 0.5 * pi + zeros(2, angle_count), a);
        else
            [r_vals, angle_vals] = bl_to_sphere(r_vals, angle_vals, a);
            angle_vals = angle_vals + 0.5 * pi;
        end
    end
    h_line_type = '-';
    if (isnan(r_vals(2)))
        h_line_type = '--';
    end
    polarplot([angle_vals(1, :), NaN, angle_vals(1, :) + pi], [r_vals(1, :), NaN, r_vals(1, :)], ['r', h_line_type]);
    polarplot([angle_vals(2, :), NaN, angle_vals(2, :) + pi], [r_vals(2, :), NaN, r_vals(2, :)], ['b', h_line_type]);
end
if (drawing_surfaces(2))
    angle_vals = ones(2, 1) * init_angle_series;
    if (is_equatorial)
        if (is_bl)
            r_vals = [0.0; 2.0 * M] * ones(1, angle_count);
        else
            r_vals = [abs(a); sqrt(4.0 * M^2.0 + a^2.0)] * ones(1, angle_count);
        end
    else
        if (abs(a) <= M)
            r_vals = calc_static_surf(angle_vals(1, :), M, a);
        else
            compl_count = 5;
            crit_theta = [acos(M / abs(a)), acos(-M / abs(a))];
            ss_i = find((init_angle_series > crit_theta(1)) & (init_angle_series < crit_theta(2)));
            angle_vals = ones(2, 1) * [linspace(crit_theta(1), init_angle_series(ss_i(1)), compl_count), ...
                                       init_angle_series(ss_i(2: 1: (end - 1))), ...
                                       linspace(init_angle_series(ss_i(end)), crit_theta(2), compl_count)];
            r_vals = [[M; M], calc_static_surf(angle_vals(1, 2: 1: (end - 1)), M, a), [M; M]];
        end
        if (~is_bl)
            [r_vals, angle_vals] = bl_to_sphere(r_vals, angle_vals, a);
        end
        angle_vals = angle_vals + 0.5 * pi;
    end
    if (solid_stat_surf)
        ss1_line_spec = 'm-'; ss2_line_spec = 'g-';
        ss_line_width = 0.5;
    else
        ss1_line_spec = 'r:'; ss2_line_spec = 'b:';
        ss_line_width = 1.0;
    end
    polarplot([angle_vals(1, :), NaN, angle_vals(1, :) + pi], [r_vals(1, :), NaN, r_vals(1, :)], ss1_line_spec, 'LineWidth', ss_line_width);
    polarplot([angle_vals(2, :), NaN, angle_vals(2, :) + pi], [r_vals(2, :), NaN, r_vals(2, :)], ss2_line_spec, 'LineWidth', ss_line_width);
end
if (drawing_surfaces(3))
    if (is_equatorial)
        if (is_bl)
            polarplot(0.0, 0.0, 'k.', 'MarkerSize', 8.0);
        else
            polarplot([init_angle_series, NaN, init_angle_series + pi], abs(a) + zeros(1, 2 * angle_count + 1), 'k-', 'LineWidth', 1.0);
        end
    else
        if (is_bl)
            polarplot(0.0, 0.0, 'k.', 'MarkerSize', 8.0);
        else
            polarplot([0.0, NaN, pi], abs(a) + zeros(1, 3), 'k.', 'MarkerSize', 8.0);
        end
    end
end
hold off;
if (~is_equatorial)
    thetaticklabels(num2str(phi_deg_to_theta_deg(str2num(char(thetaticklabels()))))); %#ok<ST2NM>
end
if (r_axis_end > 0.0)
    axis_now = axis();
    axis_now(4) = r_axis_end;
    axis(axis_now);
end
set(get(the_canvas, 'Children'), 'FontName', 'Times New Roman');

function modified_degs = phi_deg_to_theta_deg(the_degrees)
    modified_degs = 360.0 - mod(the_degrees - 90.0, 360.0);
    modified_degs(modified_degs == 360.0) = 0.0;
end
