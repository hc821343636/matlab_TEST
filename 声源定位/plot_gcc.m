clc; clear all; % �������ں����б���

sample_freq = 16000; % ���ò���Ƶ��Ϊ44100Hz
speed_sound = 343.0; % ��������Ϊ343.0m/s
K = 6; % ������˷�����Ϊ6

% ������Ƶ�ļ�·��
Audio_path = "C:\Users\82134\Desktop\";

%{
mic1=[0,0.036,0]
mic2=[-0.0312,0.018,0]
mic3=[-0.0312,-0.018,0]
mic4=[0,-0.036,0]
mic5=[0.0312,-0.018,0]
mic5=[0.0312,0.018,0]

%}
% ��ȡ������˷����Ƶ����
[T(:,1), fs] = audioread(Audio_path + "hello outside 1.wav");
[T(:,2), fs] = audioread(Audio_path + "hello outside 2.wav");
[T(:,3), fs] = audioread(Audio_path + "hello outside 3.wav");
[T(:,4), fs] = audioread(Audio_path + "hello outside 4.wav");
[T(:,5), fs] = audioread(Audio_path + "hello outside 5.wav");
[T(:,6), fs] = audioread(Audio_path + "hello outside 6.wav");

% ������˷�����
mic_coordinate(1,:) = [0 0.036 0]; % 1����˷�
mic_coordinate(2,:) = [-0.0312 0.018 0];
mic_coordinate(3,:) = [-0.0312 -0.018 0];
mic_coordinate(4,:) = [0 -0.036 0];
mic_coordinate(5,:) = [0.0312 -0.018 0];
mic_coordinate(6,:) = [0.0312 0.018 0];

% ���Ƶ�һ����Ƶͨ����ʱ������ͼ
figure;
t1 = 1:length(T(:,1));
t2 = t1 / fs;
plot(t2, T(:,1), 'LineWidth', 2, 'Color', 'b');
hold on;

% ���������ռ�ı߽�
lsb = [-2 -1 -1];
usb = [2 1 1];
%{
lsb = [-2 0 -2]��������������������ռ��һ���ǵ㣬��Ϊ��С��x��y��z����ֵ������������У���ָ����x�������СֵΪ-2�ף�y�������СֵΪ0�ף�z�������СֵΪ-2�ס�
usb = [2 2 2]��������������������ռ����һ���ǵ㣬��Ϊ����x��y��z����ֵ������������У���ָ����x��������ֵΪ2�ף�y��������ֵΪ2�ף�z��������ֵΪ2�ס�
%}
% �����һ���ŵ������ݣ�Ѱ�ҷ�ֵ
data1 = T(:,1);
ref_pk = max(data1) * 0.6;
[pks, peak_loc_total] = findpeaks(data1, 'MinPeakHeight', ref_pk, 'MinPeakDistance', 0.2 * fs);

% ��ȡ��ֵ����
[n1, n2] = size(peak_loc_total);
result = []; % ��ʼ���������
for i = 1:n1
    s = zeros(36000, K); % ��ʼ��һ�����ڴ洢��ƵƬ�εľ���
    for data_num = 1:K
        index_tmp = peak_loc_total(i); % ��ȡ��ǰ��ֵλ��
        s(:, data_num) = T(index_tmp : index_tmp + 35999, data_num); % ��ȡ��ƵƬ��
    end

    % ʹ��srppolar��srplems����������Դ��λ
    %[finalpos, finalsrp, finalfe] = srppolar(s, mic_coordinate, fs, lsb, usb); % ʹ��srppolar������λ��Դ
    [finalpos, finalsrp, finalfe] = srplems(s, mic_coordinate, fs, lsb, usb); % ʹ��srplems������λ��Դ
    x = finalpos(1); % ��ȡ��λ�õ���x����
    y = finalpos(2); % ��ȡ��λ�õ���y����
    r = sqrt(x * x + y * y); % ���㵽ԭ��ľ���
    result = [result acos(x / r) * 180 / pi]; % ���㲢�洢��λ��
    fprintf("varphi: %f\n", acos(x / r) * 180 / pi); % �����λ��

end
result