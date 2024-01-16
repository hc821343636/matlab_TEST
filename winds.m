% 数据准备
%data = % 在这里输入你的数据
function [peak_data] = winds(x,smoothData_data_medfilt,peak_threshold,y_LimitationDown)

% 初始化变量
data=smoothData_data_medfilt;
smooth_data = zeros(size(data)); % 存储光滑信号
peak_data = zeros(size(data)); % 存储波峰信号
smooth_indices = []; % 存储光滑信号的下标
peak_indices = []; % 存储波峰信号的下标

% 参数设置
window_size = 10; % 滑动窗口大小，用于平滑信号
%peak_threshold = 0.8; % 波峰阈值，用于检测波峰

% 平滑信号
for i = 1:size(data, 1)
    if i < (window_size + 1) / 2 || i > size(data, 1) - (window_size - 1) / 2
        % 边界处直接赋值
        smooth_data(i) = data(i);
    else
        % 在窗口内进行平均值计算
        smooth_data(i) = mean(data(i - round((window_size - 1) / 2):round(i + (window_size - 1) / 2)));
    end
end

% 检测波峰及其周边信号
for i = 2:size(smooth_data, 1) - 1
    if smooth_data(i) > smooth_data(i - 1) && smooth_data(i) > smooth_data(i + 1) && smooth_data(i) > mean(smooth_data) + peak_threshold * std(smooth_data)
        % 当前点为波峰
        peak_indices = [peak_indices,i];
        % 寻找波峰开始点
        start_index = i - 1;
        while start_index >= 1 && smooth_data(start_index) < smooth_data(start_index + 1)
            start_index = start_index - 1;
        end
        % 寻找波峰结束点
        end_index = i + 1;
        while end_index <= size(smooth_data, 1) && smooth_data(end_index) < smooth_data(end_index - 1)
            end_index = end_index + 1;
        end
        if start_index ==0
            start_index=1;
        end
        if end_index>size(data)
            end_index=size(data);
        end
        peak_data(start_index:end_index) = data(start_index:end_index); % 将波峰及其周边信号保存到peak_data中
    else
        % 当前点非波峰
        smooth_indices = [smooth_indices,i];
        smooth_data(i) = mean(smooth_data(i - 1:i + 1)); % 将当前点平滑处理
    end
end

% 可视化结果
figure;
hold on;
%plot(data, 'k');
%plot(smooth_data, 'b');
plot(x,peak_data, 'r');
ylim([max(peak_data)-y_LimitationDown max(peak_data)]);
legend('波峰信号');



