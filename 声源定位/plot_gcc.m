clc; clear all; % �������ں����б���

sample_freq = 16000; % ���ò���Ƶ��Ϊ44100Hz
speed_sound = 343.0; % ��������Ϊ343.0m/s
K = 6; % ������˷�����Ϊ6

% ������Ƶ�ļ�·��
Audio_path = "C:\Users\82134\Desktop\matlab_TEST-master\matlab_TEST-master\voice\6mic_10_google ";
%{  
 0-180    ����-5
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


0-360 ������
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
%ò�����������ҵߵ���
% ��ȡ������˷����Ƶ����
[T(:,6), fs] = audioread(Audio_path + "1.wav");
[T(:,5), fs] = audioread(Audio_path + "2.wav");
[T(:,4), fs] = audioread(Audio_path + "3.wav");
[T(:,3), fs] = audioread(Audio_path + "4.wav");
[T(:,2), fs] = audioread(Audio_path + "5.wav");
[T(:,1), fs] = audioread(Audio_path + "6.wav");

% ������˷�����
%{

mic_coordinate(1,:) = [0 0.036 0]; % 1����˷�
mic_coordinate(2,:) = [-0.0312 0.018 0];
mic_coordinate(3,:) = [-0.0312 -0.018 0];
mic_coordinate(4,:) = [0 -0.036 0];
mic_coordinate(5,:) = [0.0312 -0.018 0];
mic_coordinate(6,:) = [0.0312 0.018 0];
%}



mic_coordinate(6,:) = [0 0.036 0]; % 1����˷�
mic_coordinate(5,:) = [-0.036 0.0311 0];
mic_coordinate(4,:) = [-0.036 -0.0311 0];
mic_coordinate(3,:) = [0 -0.036 0];
mic_coordinate(2,:) = [0.036 -0.0311 0];
mic_coordinate(1,:) = [0.036 0.0311 0];


% ���Ƶ�һ����Ƶͨ����ʱ������ͼ
figure;
t1 = 1:length(T(:,1));
t2 = t1 / fs;
plot(t2, T(:,1), 'LineWidth', 2, 'Color', 'b');
hold on;

% ���������ռ�ı߽�
lsb = [-1 -1 0.3];
usb = [1 1 0.3];
%{
lsb = [-2 0 -2]��������������������ռ��һ���ǵ㣬��Ϊ��С��x��y��z����ֵ������������У���ָ����x�������СֵΪ-2�ף�y�������СֵΪ0�ף�z�������СֵΪ-2�ס�
usb = [2 2 2]��������������������ռ����һ���ǵ㣬��Ϊ����x��y��z����ֵ������������У���ָ����x��������ֵΪ2�ף�y��������ֵΪ2�ף�z��������ֵΪ2�ס�
%}
% �����һ���ŵ������ݣ�Ѱ�ҷ�ֵ
data1 = T(:,1);
ref_pk = max(data1) * 0.5;
[pks, peak_loc_total] = findpeaks(data1, 'MinPeakHeight', ref_pk, 'MinPeakDistance', 0.2 * fs);

% ��ȡ��ֵ����
[n1, n2] = size(peak_loc_total);
result = []; % ��ʼ���������
voice_len=32768/2;
%{
for i = 1:n1
    s = zeros(18000, K); % ��ʼ��һ�����ڴ洢��ƵƬ�εľ���
    for data_num = 1:K
        index_tmp = peak_loc_total(i); % ��ȡ��ǰ��ֵλ��
        s(:, data_num) = T(index_tmp : index_tmp + 17999, data_num); % ��ȡ��ƵƬ��
    end
 %}
for i = 1:n1
    s = zeros(voice_len,  K); % ��ʼ��һ�����ڴ洢��ƵƬ�εľ���
    for data_num = 1:K
        index_tmp = peak_loc_total(i); % ��ȡ��ǰ��ֵλ��
        start_index = max(1, index_tmp - voice_len/2+1); % ��ȡ����Ƭ�ε���ʼλ�ã���֤��Խ��
        end_index = min(length(T), index_tmp + voice_len/2); % ��ȡ����Ƭ�εĽ���λ�ã���֤��Խ��
        s(:, data_num) = T(start_index:end_index, data_num); % ��ȡ��������Ƭ��
    end
    % ʹ��srppolar��srplems����������Դ��λ
    [finalpos, finalsrp, finalfe] = srppolar(s, mic_coordinate, fs, lsb, usb); % ʹ��srppolar������λ��Դ
    [finalpos, finalsrp, finalfe] = srplems(s, mic_coordinate, fs, lsb, usb); % ʹ��srplems������λ��Դ
    x = finalpos(1); % ��ȡ��λ�õ���x����
    y = finalpos(2); % ��ȡ��λ�õ���y����
    r = sqrt(x * x + y * y); % ���㵽ԭ��ľ���

    %result = [result acos(x / r) * 180 / pi]; % ���㲢�洢��λ��  ֻ�ܻ��0-180�� 
     % ���㷽λ�ǿ��Ի��0-360��
    theta = atan2(y, x) * 180 / pi;

    % ȷ����λ����0��360�ȷ�Χ��
    if theta < 0
        theta = theta + 360;
    end
    % �洢������
    result = [result theta]; % ��������ķ�λ����ӵ�result������
    fprintf("x: %f y:%f\n", x,y);
    %fprintf("varphi: %f\n", acos(x / r) * 180 / pi); % �����λ��


end
fprintf("mean is %f",mean(result))
%result