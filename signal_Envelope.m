
n=-5000:20:5000;            % 样点设置
data=detrend(data,2,'Continuous',false);
N=length(data);                % 信号样点数
nt=0:N-1;                   % 设置样点序列号
x=120+96*exp(-(n/1500).^2).*cos(2*pi*n/600); % 设置信号
[up,down] = envelope(data,10,'peaks');
% 作图
plot(nt,data,'k',nt,up,'r',nt,down,'g');
xlabel('样点'); ylabel('幅值'); grid;
title('调用envelope函数求取上下包络曲线图')
set(gcf,'color','w');

% 生成一个示例信号
t = 0:0.1:10;
x = data;

% 使用findpeaks函数找到信号的波峰
[~, peak_locs] = findpeaks(x);

% 计算波峰的平均值
avg_peak = mean(x(peak_locs));

% 将平均值重复为与信号长度相同的向量，并将其减去信号
x_new = repmat(avg_peak, size(x)) - x;

% 绘制原始信号和处理后的信号
plot (x);
hold on;
plot( x_new)
legend('原始信号', '处理后的信号');
