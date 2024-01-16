% 以下数字除以20为秒
%405470-406470
%105470-107070
%155470-160770
%170000-171000
%275000-285000
%285000-286800
gap = 1000;%取多少秒
sample =20;%样本一秒采样20个点
row_gap=5;%一秒取5个点
y_LimitationDown=0.5;%波峰作图的显示下界

%285400 286700 286900 319600
start_n = 313000;%区间开始秒
end_n =  start_n+gap;%区间结束秒
% 提取每行第2个到最后一个元素，每
data = table2array(myfile(start_n:sample/row_gap:end_n,2:end));
peak_threshold=0.8;%波峰阈值

data_all=table2array(myfile(:, 2:1:end));

% 将数据转换为一维数组
data = reshape(data', [], 1);
x = 0:1/row_gap:(length(data)-1)/row_gap;%横坐标
data_all=reshape(data_all',[],1);

windowSize=10;%滑动窗口大小
smoothData_1 = smooth(data,windowSize);
smoothData_2=smooth(smoothData_1,windowSize);


% 设计一个低通滤波器，截止频率为200Hz
fs = 2000; % 采样率
fc = 200; % 截止频率
[b,a] = butter(6,fc/(fs/2),'low');

% 使用滤波器滤去信号的高频部分

% 绘制原始信号和滤波后的信号

t = (0:length(data)-1)/fs;
figure;
subplot(4,1,1);
data_1=data(100:end);

plot(x,data,'b','LineWidth',1);
filtered_data = filter(b,a,data);
filtered_data_exp100=filtered_data(100:end);
x_1=100/row_gap:1/row_gap:(length(filtered_data_exp100)-1+100)/row_gap;%横坐标

hold on;
plot(x,filtered_data,'r','LineWidth',1);
smoothData_3=smooth(filtered_data,windowSize);
smoothData_4=smooth(smoothData_3,windowSize);
smoothData_5=smooth(smoothData_4,windowSize);
plot(x,smoothData_5,'g','LineWidth',1)
title('滤波后的信号');
xlabel('样本');
ylabel('幅值')
legend('原始信号','低通滤波','平均滤波');
xlim([5 max(x)]);



subplot(4, 1, 2);

plot(x,data,'b','LineWidth',1);
hold on;
plot(x,smoothData_2,'r','LineWidth',2);
xlabel('秒');
ylabel('幅值');
legend('原始信号','两次平均滤波');

subplot(4, 1, 3);
plot(x,smoothData_5,'g','LineWidth',1);
hold on;
plot(x,smoothData_2,'r','LineWidth',1);

title('Original Signal');

subplot(4, 1, 4);
[b, a] = butter(7, 0.2, 'low');
data_butter = filter(b, a, data);
plot(data);
hold on;
plot(data_butter,'k');
smoothData_data_butter=smooth(data_butter,windowSize);
smoothData_data_butter=smooth(data_butter,windowSize);
hold on;
plot(smoothData_data_butter,'r');
title('Butterworth Filter');
legend('Original Signal', 'Filtered Signal');
xlim([25, 420]);


data_medfilt = medfilt1(data, 5);
figure;
plot(x,data);
hold on;
plot(x,smoothData_2,'b','LineWidth',1);
smoothData_data_medfilt=smooth(wavelet_transformation(data),windowSize);
%smoothData_data_medfilt=smoothData_data_medfilt(1:end-4);
hold on;
plot(x,smoothData_data_medfilt,'r','LineWidth',1);
legend('原始信号', '平滑滤波','小波去噪');

peak_data=winds(x,smoothData_data_medfilt,peak_threshold,y_LimitationDown);





%{
figure;
% 生成横坐标向量
%x = 0:0.05:(length(data_all)-1)*0.05;
x = 0:0.2:(length(data_all)-1)*0.2;

% 绘制折线图
plot(x, data_all);
d=datetime(1593046814, 'ConvertFrom', 'posixtime' ,'TimeZone', 'local')
%}