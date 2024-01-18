clc;clear all;
theta_ = 5:90;
varphi_ = 0:360;
theta = theta_ * pi / 180;
% theta = log(tan(pi/4+ theta_ *pi/360));
varphi = varphi_ * pi /180;
M = length(theta_);%俯仰角
N = length(varphi_);%水平偏向角
u = sin(theta)'*cos(varphi);
v = sin(theta)'*sin(varphi);
w = repmat(cos(theta)',1,N);
sample_freq = 44100;
speed_sound = 343.0;
Audio_path = "D:\matlab\matlabR2019b\bin\microphone_array\sound_test_4\office_44_1K_sampling_2-3K_Fre_5ms_duration_5s_Inter\1\";




%麦克风标号与读取的wav文件标号似乎是上下左右颠倒的
[T(:,1),fs] = audioread(Audio_path+"Audio Track-4.wav");
[T(:,2),fs] = audioread(Audio_path+"Audio Track-5.wav");
[T(:,3),fs] = audioread(Audio_path+"Audio Track-6.wav");
[T(:,4),fs] = audioread(Audio_path+"Audio Track.wav");
[T(:,5),fs] = audioread(Audio_path+"Audio Track-2.wav");
[T(:,6),fs] = audioread(Audio_path+"Audio Track-3.wav");

% mic_coordinate(1,:) = [0      0   0]; %1# microphone
% mic_coordinate(2,:) = [0.05   0   0];
% mic_coordinate(3,:) = [0.075  0   0.0425];
% mic_coordinate(4,:) = [0.05   0   0.085];
% mic_coordinate(5,:) = [0      0   0.085];
% mic_coordinate(6,:) = [-0.025 0   0.0425];   
%水平放置的坐标系
mic_coordinate(1,:) = [0      0   0]; %1# microphone
mic_coordinate(2,:) = [0.05   0   0];
mic_coordinate(3,:) = [0.075  -0.0425  0];
mic_coordinate(4,:) = [0.05   -0.085   0];
mic_coordinate(5,:) = [0      -0.085   0];
mic_coordinate(6,:) = [-0.025 -0.0425  0]; 

search_coordinate(:,:,1)=u;
search_coordinate(:,:,2)=v;
search_coordinate(:,:,3)=w;

%%%Parameters for chirp signals
chirp_length = 0.002;  %%%the length of chirp 2 ms
chirp_interval = 0.5;  %%%the interval of chirp
chirp_freq1 = 2000; 
chirp_freq2 = 3000; 
win_beam = round((0.012)*fs);
start_beam = round(0.002*fs);

K = size(mic_coordinate,1);
angle_cos = zeros(M,N,K);

% figure
% t1 = 1:length( T(:,1) );
% t2 =  t1/fs;
% plot(t2,  T(:,1),'LineWidth',2,'Color','b');
% set(gca,'ytick',[]) 
% set(gca,'xtick',[]) 
% hold on;

tt = 0:1/fs:chirp_length;
sampleSignal = chirp(tt,chirp_freq1,chirp_length,chirp_freq2,'linear');
T = bandpass(T,[2000 3000],fs);

% figure
% t1 = 1:length( T(:,1) );
% t2 =  t1/fs;
% plot(t2,  T(:,1),'LineWidth',2,'Color','b');
% set(gca,'ytick',[]) 
% set(gca,'xtick',[]) 
% hold on;


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
Energy_matrix2 = zeros(M,N);
%为了获得len_E的笨办法
 waves = {};
 index_tmp = peak_loc_total(2, 1);
 for data_num = 1:n2
    data_seg = T(index_tmp-round(0.008*fs): index_tmp+win_beam,data_num);
    waves{data_num} = data_seg'; 
 end
len_E = length(delay_and_sum_beamforming_v5_yanzhi(waves,delay_point,start_beam,M,N,1));
Energy_matrix = zeros(M,N,len_E);
count = 0;

cnt = 1;
distance = 0; distance1 = 0;
theta_ave = 0;varphi_ave = 0;
true_cnt = 0;
save_path = 'D:\microphone\dataset\dataset1\image1.';
parameter = 2.3;
for cnt = 2:4:20

    for chirp_num = cnt:cnt+3
%         waves = {};
        index_tmp = peak_loc_total(chirp_num, 1);
        for data_num = 1:n2
            %用于Energy_matrix,包含了直达信号
        data_seg1 = T(index_tmp-round(0.008*fs): index_tmp+win_beam,data_num);
        waves1{data_num} = data_seg1'; 
        %用于Energy_matrix2，不包含直达信号
        data_seg2 = T(index_tmp: index_tmp+win_beam,data_num);
        waves2{data_num} = data_seg2'; 
        end
    
        %Energy_matrix = delay_and_sum_beamforming_v5_yanzhi(waves,delay_point,start_beam,M,N,K);
        Energy_matrix = Energy_matrix + delay_and_sum_beamforming_v5_yanzhi(waves1,delay_point,start_beam,M,N,K);
        Energy_matrix2 = Energy_matrix2 + delay_and_sum_beamformingv4(waves2,delay_point,start_beam,M,N,K);
    end
    Energy_matrix = Energy_matrix/4;
    Energy_matrix2 = Energy_matrix2/4;

%     surf(u,v,w,Energy_matrix2, 'edgecolor', 'none');
%     set(gca,'ytick',[]) 
%     set(gca,'xtick',[])  
%     set(gca,'ztick',[])
%     hidden off;

    %进行墨卡托投影
    M_mokato = ceil(log(tan(pi/4+ 85 *pi/360)) * 180 / pi);
    Energy_matrix_mokato = zeros(M_mokato+ M ,N);
    row = 1;
    for i = theta_
        if i + 1 > 85
            break;
        end
        diff = ceil( (log(tan(pi/2 - i * pi/360))-log(tan(pi/2 - (i+1) * pi/360))) * 180 / pi );
        while diff > 0
            Energy_matrix_mokato(row,:) = Energy_matrix2(i, :);
            diff = diff - 1;
            row = row + 1;
        end
    end
    Energy_matrix_mokato = Energy_matrix_mokato(1:row-1 , :);
    figure('visible','off');
%     figure;
    image(Energy_matrix_mokato,'CDataMapping','scaled');
    colorbar('off');
%    caxis([-1,1]);
    set(gca,'ytick',[])  %隐去y轴坐标值
    set(gca,'xtick',[])  %隐去x轴坐标值
    set(gca,'Position',[ 0 0 1 1]);
    box off;
    grid on;

    name_string = [ save_path  num2str(count)  '.png'];
    saveas(gcf,name_string);
    count = count + 1;
    close all;
    
    
    %data augmentation
%     for k = 1.2:0.2:1.7
%         E_attenuate = real(10.^(0.8*(10 * log10(Energy_matrix2) - 20*log10(k))/10));
%         figure('visible','off');
%         image(E_attenuate,'CDataMapping','scaled');
%         colorbar('off');
%         caxis([min_val,max_val]);
%     set(gca,'ytick',[])  %隐去y轴坐标值
%     set(gca,'xtick',[])  %隐去x轴坐标值
%     box off;
%     grid on;
%     name_string = [ save_path  num2str(count)  '.png'];
%     saveas(gcf,name_string);
%     count = count + 1;
%     end
    

end
