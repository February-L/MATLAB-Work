clear;
% ����Bardeenʱ����Ч�ƣ������theta = pi/2��
% ��Ȼ��λ�ƣ�c = G = 1, mu_0 = 4*pi

global M; % ����Դ��������ʷ����ʱ�ն����ӽ�λ��r = 2M��
global g; % ���򻯲������ɽ���Ϊһ�����Դŵ����ӵĴźɣ���(g/M)^2 < 16/27ʱBardeenʱ�մ��������ӽ�
global const_L; % �������ת����phi���غ�����L = r^2 * d(phi)/ds��������˶�s��Ϊ����������
global motion_timelike; motion_timelike = false; % �Ƿ�Ϊ��ʱ�˶���������Ϊ����˶�
M = 1.0;
g = 0.7; const_L = 3.6;
var_choice = 'L'; % 'L' - �仯L�̶�g��'g' - �仯g�̶�L
%var_vals = [0.0, 0.7, 0.85, 0.9, 1.0];
%var_vals = [3.0, 3.2, 3.4, 3.6];
val_vals = [];
var_count = length(var_vals); % ��var_countΪ0�򲻶�g��L��var_vals���¸�ֵ\

r_lower = 0.0; r_upper = 16.0;
r_ticks = r_lower: 0.01: r_upper;
r_vals = r_ticks;
%r_ticks = sqrt(r_lower): 0.01: sqrt(r_upper);
%r_vals = r_ticks .^ 2.0;
V2_lower = 0.0; V2_upper = 1.0; % V^2�������������ƣ���ֵΪ[]�򲻲�ȡ��Ӧ����

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
            case 'g'
                g = var_vals(i);
            otherwise
                fprintf("Error: unavailable variable choice.\n");
                return;
        end
        var_names{i} = sprintf("\\fontsize{10}\\it%c\\rm = %g", var_choice, var_vals(i));
        plot(r_ticks, effp_square(r_vals));
    else
        plot(r_ticks, effp_square(r_vals), "k-");
    end
end
%plot(r_ticks, zeros(1, length(r_ticks)), "k:");
axis_now = axis();
if (~isempty(V2_lower))
    axis_now(3) = V2_lower;
end
if (~isempty(V2_upper))
    axis_now(4) = V2_upper;
end
axis(axis_now);
if (var_count ~= 0)
    legend(var_names);
end
xlabel("\it{r}\rm{ / }\it{M}"); ylabel("\it{V}\fontsize{4}{ }\fontsize{11}\rm{^2}", 'Rotation', 0.0);
%ylabel("-\it{g}\rm_{00}", 'Rotation', 0.0);
%xticklabels(num2str(str2num(char(xticklabels())).^2.0));
hold off;
set(get(the_canvas, 'Children'), 'FontName', 'Times New Roman');

function return_val = effp_square(r)
    % ��Ч�Ƶ�ƽ��
    global motion_timelike M g const_L;
    if (motion_timelike) % ��ʱ�˶�
        return_val = (1.0 - 2.0 .* M .* r.^2.0 ./ (r.^2.0 + g.^2.0).^1.5) .* (1.0 + (const_L./r).^2.0);
        %return_val = (1.0 - 2.0 .* M .* r.^2.0 ./ (r.^2.0 + g.^2.0).^1.5); % ��Ϊ���ضȹ�(-g00)�༴�ȹ�(1/g11)��ֵ
    else % ����˶�
        return_val = (1.0 - 2.0 .* M .* r.^2.0 ./ (r.^2.0 + g.^2.0).^1.5) .* (const_L./r).^2.0;
    end
end
