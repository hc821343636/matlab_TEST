clc;clear all;
theta_ = 10:170;
varphi_ = 10:170;
theta = theta_ * pi / 180;
varphi = varphi_ * pi /180;
M = length(theta);%俯仰角
N = length(varphi);%水平偏向角
u = sin(theta)'*cos(varphi);
v = sin(theta)'*sin(varphi);
w = repmat(cos(theta)',1,N);
sample_freq = 44100;
Point = sample_freq * 0.6;
speed_sound = 343.0;
Audio_path = "D:\matlab\matlabR2019b\bin\microphone_array\sound_test_4\office_44_1K_sampling_2-3K_Fre_5ms_duration_5s_Inter\1\";



% [T(:,1),fs] = audioread(Audio_path+"Audio Track-4.wav",[1000 Point]);
% [T(:,2),fs] = audioread(Audio_path+"Audio Track-5.wav",[1000 Point]);
% [T(:,3),fs] = audioread(Audio_path+"Audio Track-6.wav",[1000 Point]);
% [T(:,4),fs] = audioread(Audio_path+"Audio Track.wav",[1000 Point]);
% [T(:,5),fs] = audioread(Audio_path+"Audio Track-2.wav",[1000 Point]);
% [T(:,6),fs] = audioread(Audio_path+"Audio Track-3.wav",[1000 Point]);
%麦克风标号与读取的wav文件标号似乎是上下左右颠倒的
[T(:,1),fs] = audioread(Audio_path+"Audio Track-4.wav");
[T(:,2),fs] = audioread(Audio_path+"Audio Track-5.wav");
[T(:,3),fs] = audioread(Audio_path+"Audio Track-6.wav");
[T(:,4),fs] = audioread(Audio_path+"Audio Track.wav");
[T(:,5),fs] = audioread(Audio_path+"Audio Track-2.wav");
[T(:,6),fs] = audioread(Audio_path+"Audio Track-3.wav");

mic_coordinate(1,:) = [0      0   0]; %1# microphone
mic_coordinate(2,:) = [0.05   0   0];
mic_coordinate(3,:) = [0.075  0   0.0425];
mic_coordinate(4,:) = [0.05   0   0.085];
mic_coordinate(5,:) = [0      0   0.085];
mic_coordinate(6,:) = [-0.025 0   0.0425];   

search_coordinate(:,:,1)=u;
search_coordinate(:,:,2)=v;
search_coordinate(:,:,3)=w;

%%%Parameters for chirp signals
chirp_length = 0.002;  %%%the length of chirp 2 ms
chirp_interval = 0.5;  %%%the interval of chirp
chirp_freq1 = 2000; 
chirp_freq2 = 3000; 
win_beam = round((0.04)*fs);
start_beam = round(0.0034*fs);

K = size(mic_coordinate,1);
angle_cos = zeros(M,N,K);
figure
t1 = 1:length( T(:,1) );
t2 =  t1/fs;
plot(t2,  T(:,1),'LineWidth',2,'Color','b');
hold on;
% a = 5*fs;b =8*fs;
% num = 3;
% sm = 0;
% for i =1:6
%     sm  = sm + sum(abs(T(a:b,i)));
% end
% ave_6 = sm /num/6;
% ave = sum(abs(T(a:b,1)))/num;
% fprintf("ave: %d,ave_s: %d \n",ave,ave/fs);
% fprintf("ave_6 :%d,ave_6_s:%d \n",ave_6,ave_6/fs);


% Calculate delay points
for i=2:K
    for j=1:M
        for k=1:N
            angle_cos(j,k,i)=squeeze(mic_coordinate(i,:))*squeeze(search_coordinate(j,k,:));
        end
    end
end
delay_point = round(angle_cos / speed_sound * fs);

peak_loc_total = [];
for i = 1:K
   data1 = T(:,i);
   [pks, peak_loc] = findpeaks(data1,'MinPeakHeight',0.05,'MinPeakDistance',(chirp_interval -0.02)*fs);
   peak_loc_total = [peak_loc_total peak_loc]; 
end


[n1, n2] = size(peak_loc_total); 
count = 0;
for chirp_num = 4:4
    waves = {};
    index_tmp = peak_loc_total(chirp_num, 1);
    for data_num = 1:n2
       data_seg = T(index_tmp: index_tmp+win_beam,data_num);
       waves{data_num} = data_seg'; 
    end
    
    while(start_beam<=round(0.008*fs))
        Energy_matrix = delay_and_sum_beamformingv4(waves,delay_point,start_beam,M,N,K);
    [x,y]=find(Energy_matrix==max(max(Energy_matrix)));
    theta_max = theta_(x);
    varphi_max = varphi_(y);
    fprintf("varphi: %d, theta: %d \n",theta_max,varphi_max);
    fprintf("*****************************\n");
    
    figure('visible','off');
    %figure;
    image(Energy_matrix,'CDataMapping','scaled');
    colorbar('off');
    caxis([0,0.15]);
    set(gca,'ytick',[])  %隐去y轴坐标值
    set(gca,'xtick',[])  %隐去x轴坐标值
%     xlabel('phiAngles','FontWeight', 'bold','FontSize', 15);
%     ylabel('thetaAngles','FontWeight', 'bold','FontSize', 15);
%     set(gca,'FontWeight', 'bold');
%     set(gca,'FontSize', 15);
    box off;
    grid on;
    name_string = [ 'D:\microphone\dataset\dataset2\image1.'  num2str(count)  '.png'];
    saveas(gcf,name_string);
    start_beam = start_beam + round(0.0002*fs);
    count = count + 1;
    end
   start_beam = round(0.0034*fs);
end
