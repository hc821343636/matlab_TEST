clc; clear all; % 清除命令窗口和所有变量

sample_freq = 16000; % 设置采样频率为44100Hz
speed_sound = 343.0; % 设置声速为343.0m/s
K = 6; % 设置麦克风数量为6

% 设置音频文件路径
Audio_path = "C:\Users\82134\Desktop\matlab_TEST-master\matlab_TEST-master\voice\6mic_10_google ";
%{  
 0-180    修正-5
        true    test    diff
0.5      60       63      3.
1        30        36     6.
1.5      0         0.8    0.8.
2        30        28       2.
2.5      60         61      1.
3       90          93      3.
3.5     120         130     10.
4       150         156     6.
4.5     180         178     2.
5       150          152    2.
5.5     120         135     15.
6       90           102    12.
mean                        5.2.


0-360 无修正
        true    test    diff
0.5      60       71      11
1        30        43     13
1.5      0         6      6
2        330       331    1
2.5      300        294   6
3       270          260  10
3.5     240         226  14
4       210         199  11
4.5     180         180  0
5       150          165 15
5.5     120         147  27
6       90           108 18
mean                     11.0
%}
%{
mic1=[0,0.036,0]
mic2=[-0.0312,0.018,0]
mic3=[-0.0312,-0.018,0]
mic4=[0,-0.036,0]
mic5=[0.0312,-0.018,0]
mic5=[0.0312,0.018,0]

%}
%貌似是上下左右颠倒的
% 读取六个麦克风的音频数据
[T(:,6), fs] = audioread(Audio_path + "1.wav");
[T(:,5), fs] = audioread(Audio_path + "2.wav");
[T(:,4), fs] = audioread(Audio_path + "3.wav");
[T(:,3), fs] = audioread(Audio_path + "4.wav");
[T(:,2), fs] = audioread(Audio_path + "5.wav");
[T(:,1), fs] = audioread(Audio_path + "6.wav");

% 定义麦克风坐标
%{

mic_coordinate(1,:) = [0 0.036 0]; % 1号麦克风
mic_coordinate(2,:) = [-0.0312 0.018 0];
mic_coordinate(3,:) = [-0.0312 -0.018 0];
mic_coordinate(4,:) = [0 -0.036 0];
mic_coordinate(5,:) = [0.0312 -0.018 0];
mic_coordinate(6,:) = [0.0312 0.018 0];
%}



mic_coordinate(6,:) = [0 0.036 0]; % 1号麦克风
mic_coordinate(5,:) = [-0.036 0.0311 0];
mic_coordinate(4,:) = [-0.036 -0.0311 0];
mic_coordinate(3,:) = [0 -0.036 0];
mic_coordinate(2,:) = [0.036 -0.0311 0];
mic_coordinate(1,:) = [0.036 0.0311 0];


% 绘制第一个音频通道的时间序列图
figure;
t1 = 1:length(T(:,1));
t2 = t1 / fs;
plot(t2, T(:,1), 'LineWidth', 2, 'Color', 'b');
hold on;

% 定义搜索空间的边界
lsb = [-1 -1 0.3];
usb = [1 1 0.3];
%{
lsb = [-2 0 -2]：这个向量定义了搜索空间的一个角点，作为最小的x、y和z坐标值。在这个例子中，它指定了x坐标的最小值为-2米，y坐标的最小值为0米，z坐标的最小值为-2米。
usb = [2 2 2]：这个向量定义了搜索空间的另一个角点，作为最大的x、y和z坐标值。在这个例子中，它指定了x坐标的最大值为2米，y坐标的最大值为2米，z坐标的最大值为2米。
%}
% 处理第一个信道的数据，寻找峰值
data1 = T(:,1);
ref_pk = max(data1) * 0.5;
[pks, peak_loc_total] = findpeaks(data1, 'MinPeakHeight', ref_pk, 'MinPeakDistance', 0.2 * fs);

% 获取峰值数量
[n1, n2] = size(peak_loc_total);
result = []; % 初始化结果数组
voice_len=32768/2;
%{
for i = 1:n1
    s = zeros(18000, K); % 初始化一个用于存储音频片段的矩阵
    for data_num = 1:K
        index_tmp = peak_loc_total(i); % 获取当前峰值位置
        s(:, data_num) = T(index_tmp : index_tmp + 17999, data_num); % 截取音频片段
    end
 %}
for i = 1:n1
    s = zeros(voice_len,  K); % 初始化一个用于存储音频片段的矩阵
    for data_num = 1:K
        index_tmp = peak_loc_total(i); % 获取当前峰值位置
        start_index = max(1, index_tmp - voice_len/2+1); % 截取语音片段的起始位置，保证不越界
        end_index = min(length(T), index_tmp + voice_len/2); % 截取语音片段的结束位置，保证不越界
        s(:, data_num) = T(start_index:end_index, data_num); % 截取整个语音片段
    end
    % 使用srppolar和srplems函数进行声源定位
    [finalpos, finalsrp, finalfe] = srppolar(s, mic_coordinate, fs, lsb, usb); % 使用srppolar函数定位声源
    [finalpos, finalsrp, finalfe] = srplems(s, mic_coordinate, fs, lsb, usb); % 使用srplems函数定位声源
    x = finalpos(1); % 获取定位得到的x坐标
    y = finalpos(2); % 获取定位得到的y坐标
    r = sqrt(x * x + y * y); % 计算到原点的距离

    %result = [result acos(x / r) * 180 / pi]; % 计算并存储方位角  只能获得0-180度 
     % 计算方位角可以获得0-360度
    theta = atan2(y, x) * 180 / pi;

    % 确保方位角在0到360度范围内
    if theta < 0
        theta = theta + 360;
    end
    % 存储计算结果
    result = [result theta]; % 将计算出的方位角添加到result数组中
    fprintf("x: %f y:%f\n", x,y);
    %fprintf("varphi: %f\n", acos(x / r) * 180 / pi); % 输出方位角


end
fprintf("mean is %f",mean(result))
%result