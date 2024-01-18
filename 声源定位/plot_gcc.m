clc; clear all; % 清除命令窗口和所有变量

sample_freq = 16000; % 设置采样频率为44100Hz
speed_sound = 343.0; % 设置声速为343.0m/s
K = 6; % 设置麦克风数量为6

% 设置音频文件路径
Audio_path = "C:\Users\82134\Desktop\";

%{
mic1=[0,0.036,0]
mic2=[-0.0312,0.018,0]
mic3=[-0.0312,-0.018,0]
mic4=[0,-0.036,0]
mic5=[0.0312,-0.018,0]
mic5=[0.0312,0.018,0]

%}
% 读取六个麦克风的音频数据
[T(:,1), fs] = audioread(Audio_path + "hello outside 1.wav");
[T(:,2), fs] = audioread(Audio_path + "hello outside 2.wav");
[T(:,3), fs] = audioread(Audio_path + "hello outside 3.wav");
[T(:,4), fs] = audioread(Audio_path + "hello outside 4.wav");
[T(:,5), fs] = audioread(Audio_path + "hello outside 5.wav");
[T(:,6), fs] = audioread(Audio_path + "hello outside 6.wav");

% 定义麦克风坐标
mic_coordinate(1,:) = [0 0.036 0]; % 1号麦克风
mic_coordinate(2,:) = [-0.0312 0.018 0];
mic_coordinate(3,:) = [-0.0312 -0.018 0];
mic_coordinate(4,:) = [0 -0.036 0];
mic_coordinate(5,:) = [0.0312 -0.018 0];
mic_coordinate(6,:) = [0.0312 0.018 0];

% 绘制第一个音频通道的时间序列图
figure;
t1 = 1:length(T(:,1));
t2 = t1 / fs;
plot(t2, T(:,1), 'LineWidth', 2, 'Color', 'b');
hold on;

% 定义搜索空间的边界
lsb = [-2 -1 -1];
usb = [2 1 1];
%{
lsb = [-2 0 -2]：这个向量定义了搜索空间的一个角点，作为最小的x、y和z坐标值。在这个例子中，它指定了x坐标的最小值为-2米，y坐标的最小值为0米，z坐标的最小值为-2米。
usb = [2 2 2]：这个向量定义了搜索空间的另一个角点，作为最大的x、y和z坐标值。在这个例子中，它指定了x坐标的最大值为2米，y坐标的最大值为2米，z坐标的最大值为2米。
%}
% 处理第一个信道的数据，寻找峰值
data1 = T(:,1);
ref_pk = max(data1) * 0.6;
[pks, peak_loc_total] = findpeaks(data1, 'MinPeakHeight', ref_pk, 'MinPeakDistance', 0.2 * fs);

% 获取峰值数量
[n1, n2] = size(peak_loc_total);
result = []; % 初始化结果数组
for i = 1:n1
    s = zeros(36000, K); % 初始化一个用于存储音频片段的矩阵
    for data_num = 1:K
        index_tmp = peak_loc_total(i); % 获取当前峰值位置
        s(:, data_num) = T(index_tmp : index_tmp + 35999, data_num); % 截取音频片段
    end

    % 使用srppolar和srplems函数进行声源定位
    %[finalpos, finalsrp, finalfe] = srppolar(s, mic_coordinate, fs, lsb, usb); % 使用srppolar函数定位声源
    [finalpos, finalsrp, finalfe] = srplems(s, mic_coordinate, fs, lsb, usb); % 使用srplems函数定位声源
    x = finalpos(1); % 获取定位得到的x坐标
    y = finalpos(2); % 获取定位得到的y坐标
    r = sqrt(x * x + y * y); % 计算到原点的距离
    result = [result acos(x / r) * 180 / pi]; % 计算并存储方位角
    fprintf("varphi: %f\n", acos(x / r) * 180 / pi); % 输出方位角

end
result