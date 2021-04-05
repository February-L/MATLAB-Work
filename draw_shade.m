clear;
% ª≠“ªøÈ“ı”∞

shade_step = 0.1;
the_func = @(x, x0) (x - x0);
x_vals = -10.0: shade_step: 10.0;
x0_vals = [0.0, -10.0];
line_count = length(x_vals);
line_rgb = [0.5, 0.5, 0.5];

the_canvas = figure();
hold on;
for i = 1: 1: line_count
    plot(the_func(x_vals(i), x0_vals), [0.0, 10.0], 'Color', line_rgb);
end
hold off;
xlim([0.0, 10.0]);
set(get(the_canvas, 'Children'), 'Box', false, 'XColor', [1.0, 1.0, 1.0], 'YColor', [1.0, 1.0, 1.0]);
xticklabels(strings(1, length(xticklabels())));
yticklabels(strings(1, length(yticklabels())));
