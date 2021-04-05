clear;
% 绘制Kerr时空外部赤道面（B-L坐标theta = pi/2）度规
% 自然单位制：c = G = 1

M = 1.0; % 引力源质量
a = 0.8; % 转动参数
metric_funcs = {@(r) (-(1.0 - 2.0*M ./ r)); % g_tt
                @(r) (r.^2.0 + a^2.0 + 2.0*M*a^2.0 ./ r); % g_phiphi
                @(r) (-2.0*M*a ./ r); % g_tphi
                @(r) (r.^2.0 ./ (r.^2.0 - 2.0*M .* r + a^2.0)); % g_rr
                @(r) (r.^2.0 - 2.0*M .* r + a^2.0)}; % horizon_omega
drawing_curves = logical([1; 1; 1; 1; 0]); % 决定是否绘制对应度规/函数
metric_names = ["\it{g_{tt}}"; "\it{g_{\phi\phi}}"; "\it{g_{t\phi}}"; "\it{g_{rr}}"; "\it\Omega"];
% 第1~5个元素分别对应g_tt、g_phiphi、g_tphi、g_rr、horizon_omega

r_lower = 0.0; r_upper = 5.0;
y_lower = -10.0; y_upper = 10.0; % y坐标轴上下限制，若值为空则不采取对应限制

the_canvas = figure();
plot([r_lower, r_upper], zeros(1, 2), 'k-');
hold on;
func_lines = fplot(metric_funcs(drawing_curves), [r_lower, r_upper]);
legend(func_lines, metric_names(drawing_curves));

axis_now = axis();
if (~isempty(y_lower))
    axis_now(3) = y_lower;
end
if (~isempty(y_upper))
    axis_now(4) = y_upper;
end
axis(axis_now);
xlabel("\it{r}"); ylabel("\it{g_{\mu\nu}}", 'Rotation', 0.0);
set(get(the_canvas, 'Children'), 'FontName', 'Times New Roman');
set(get(the_canvas, 'Children'), 'Box', true);
