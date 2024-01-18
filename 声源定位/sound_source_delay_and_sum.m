clc;clear all;
%{
这个MATLAB脚本实现了一个基于延迟和求和（Delay-and-Sum）的波束形成算法，用于定位声源。它首先定义了一系列参数和变量，然后执行声源定位的计算。以下是对主要部分的详细解释：

初始化和参数定义:

theta_ 和 varphi_：定义搜索空间的俯仰角和水平偏向角的范围（以度为单位），然后转换为弧度。
M 和 N：计算俯仰角和水平偏向角的长度。
u, v, w：计算搜索空间中每个点的方向余弦。
sample_freq：定义采样频率。
Audio_path：定义音频文件的路径。
audioread：读取音频文件。
麦克风坐标:

mic_coordinate：定义了麦克风阵列的坐标。
搜索坐标:

search_coordinate：根据u, v, w计算搜索空间的坐标。
声波参数:

定义了chirp信号的参数，如chirp_length，chirp_interval，chirp_freq1 和 chirp_freq2。
波束形成参数和初始化:

angle_cos：预先分配用于存储角度余弦的矩阵。
绘制时间序列图:

使用plot函数绘制第一个音频通道的时间序列。
计算延迟点:

计算每个麦克风相对于每个搜索方向的延迟点。
数据处理和波束形成:

选择参考信道的峰值作为同步点。
对于每个峰值，提取各个麦克风的数据段并执行延迟和求和波束形成。
Energy_matrix：存储了每个方向上的能量值。
结果可视化:

使用image和surf函数绘制能量矩阵和三维声场分布。
声源定位:

根据能量矩阵找到最大值对应的角度，从而定位声源。
整个脚本的目的是利用多个麦克风记录的音频数据来定位声源。它通过计算不同方向上的声音能量分布，找出声源最可能的位置。
%}
theta_ = 0:180;
varphi_ = 0:180;
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
Audio_path = "C:\Users\82134\Desktop\";

%麦克风标号与读取的wav文件标号似乎是上下左右颠倒的
[T(:,1),fs] = audioread(Audio_path+"0 60 0 hello 1.wav");
[T(:,2),fs] = audioread(Audio_path+"0 60 0 hello 2.wav");
[T(:,3),fs] = audioread(Audio_path+"0 60 0 hello 3.wav");
[T(:,4),fs] = audioread(Audio_path+"0 60 0 hello 4.wav");
[T(:,5),fs] = audioread(Audio_path+"0 60 0 hello 5.wav");
[T(:,6),fs] = audioread(Audio_path+"0 60 0 hello 6.wav");


%{
mic1=[0,0.036,0]
mic2=[-0.0312,0.018,0]
mic3=[-0.0312,-0.018,0]
mic4=[0,-0.036,0]
mic5=[0.0312,-0.018,0]
mic5=[0.0312,0.018,0]

%}
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
chirp_interval = 0.7;  %%%the interval of chirp
chirp_freq1 = 2000; 
chirp_freq2 = 3000; 
win_beam = round((0.04)*fs);
start_beam = round(0.001*fs);

K = size(mic_coordinate,1);
angle_cos = zeros(M,N,K);

figure
t1 = 1:length( T(:,1) );
t2 =  t1/fs;
plot(t2,  T(:,1),'LineWidth',2,'Color','b');
hold on;

% Calculate delay points
for i=2:K
    for j=1:M
        for k=1:N
            angle_cos(j,k,i)=squeeze(mic_coordinate(i,:))*squeeze(search_coordinate(j,k,:));
        end
    end
end
delay_point = round(angle_cos / speed_sound * fs);

data1 = T(:,1); % 第一个信道的数据
[ref_pk,ref_pos] = max(data1);
[pks, peak_loc_total] = findpeaks(data1,'MinPeakHeight',ref_pk*0.7,'MinPeakDistance',(chirp_interval -0.02)*fs);

% peak_loc_total = [peak_loc_total; ref_pos];

[n1, n2] = size(peak_loc_total); 
for beep_num = 1:n1
   waves = {};
    for data_num = 1:K
%         index_tmp = peak_loc_total(beep_num,1);%应该选参考信道的peak，而不是每一个信道的peak
        index_tmp = ref_pos; % 只去能量最强的一段信号
        data_seg = T(index_tmp: index_tmp+fs,data_num);
        waves{data_num} = data_seg'; 
    end
    
    Energy_matrix = delay_and_sum_beamformingv4(waves,delay_point,start_beam,M,N,K);
    [x,y]=find(Energy_matrix==max(max(Energy_matrix)));
    theta_max = theta_(x);
    varphi_max = varphi_(y);

    fprintf("varphi: %d, theta: %d \n",varphi_max,theta_max);
    fprintf("*****************************\n");

    figure;
    image(Energy_matrix,'CDataMapping','scaled');
    colorbar('off');
    %    caxis([-1,1]);
    set(gca,'ytick',[])  %隐去y轴坐标值
    set(gca,'xtick',[])  %隐去x轴坐标值
    set(gca,'Position',[ 0 0 1 1]);
    box off;
    grid on; 
    
    figure; %三维图
    surf(u,v,w,Energy_matrix, 'edgecolor', 'none');
end
