clc;clear all;
sample_freq = 44100;
speed_sound = 343.0;
K=6;
Audio_path = "C:\Users\82134\Desktop\";

[T(:,1), fs] = audioread(Audio_path + "hello outside 1.wav");
[T(:,2), fs] = audioread(Audio_path + "hello outside 2.wav");
[T(:,3), fs] = audioread(Audio_path + "hello outside 3.wav");
[T(:,4), fs] = audioread(Audio_path + "hello outside 4.wav");
[T(:,5), fs] = audioread(Audio_path + "hello outside 5.wav");
[T(:,6), fs] = audioread(Audio_path + "hello outside 6.wav");


mic_coordinate(1,:) = [0 0.036 0]; % 1号麦克风
mic_coordinate(2,:) = [-0.0312 0.018 0];
mic_coordinate(3,:) = [-0.0312 -0.018 0];
mic_coordinate(4,:) = [0 -0.036 0];
mic_coordinate(5,:) = [0.0312 -0.018 0];
mic_coordinate(6,:) = [0.0312 0.018 0];

% mic_coordinate(1,:) = [0      0   0]; %1# microphone
% mic_coordinate(2,:) = [0.05   0   0];
% mic_coordinate(3,:) = [0.075  -0.0425  0];
% mic_coordinate(4,:) = [0.05   -0.085   0];
% mic_coordinate(5,:) = [0      -0.085   0];
% mic_coordinate(6,:) = [-0.025 -0.0425  0]; 

figure
t1 = 1:length( T(:,1) );
t2 =  t1/fs;
plot(t2,  T(:,1),'LineWidth',2,'Color','b');
hold on;


lsb=[-2 -2 -1];
usb=[2 2 1];

data1 = T(:,1); % 第一个信道的数据
ref_pk = max(data1) * 0.7;
[pks, peak_loc_total] = findpeaks(data1,'MinPeakHeight',ref_pk,'MinPeakDistance',0.2*fs);

[n1, n2] = size(peak_loc_total); 
result = [];
for i=1:n1
    s = zeros(32768,K);
    for data_num = 1:K
        index_tmp = peak_loc_total(i);
%         index_tmp = round(2.93*fs);
        s(:,data_num) = T(index_tmp:index_tmp+32767,data_num);
    end

    [finalpos,finalsrp,finalfe]=srplems(s, mic_coordinate, fs, lsb, usb);
    x=finalpos(1);
    y=finalpos(2);
    r = sqrt(x*x+y*y);
    result = [result acos(x/r)*180/pi];
    fprintf("varphi: %f\n",acos(x/r)*180/pi);
end
result
