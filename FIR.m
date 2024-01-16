% 生成随机信号
data = randn(1, 1000);
data = reshape(data', [], 1);
% 设计30点低通FIR滤波器，截止频率为0.2
b = fir1(60, 0.2);

% 应用滤波器
filtered_data = filter(b, 1, data);

% 绘制原始信号和滤波后的信号对比图
figure;
subplot(2,1,1);
plot(data);
title('原始信号');
subplot(2,1,2);
plot(filtered_data);
title('滤波后信号');
