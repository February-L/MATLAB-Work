clear;
% 以B-L坐标或球坐标绘制Kerr时空赤道面或沿转动轴截面上的视界、静止界限及奇环曲线
% 自然单位制：c = G = 1

M = 1.0; % 引力源质量
a = 1.2; % 转动参数
drawing_surfaces = logical([1, 1, 1]); % 决定是否绘制对应界面，第一、二、三个元素分别对应视界、静止界限、奇环
is_equatorial = false; % 所取绘制截面是否为赤道面，否则为沿转动轴的截面
is_bl = false; % 所取绘制坐标是否为B-L坐标，否则为球坐标
solid_stat_surf = true; % 是否将静止界限画为实线，否则画为虚点线

r_axis_end = 2.5; % 指定图像展示到r为多少（不大于0时无视此参数）
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
