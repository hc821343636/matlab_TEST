clc;clear all;
%{
���MATLAB�ű�ʵ����һ�������ӳٺ���ͣ�Delay-and-Sum���Ĳ����γ��㷨�����ڶ�λ��Դ�������ȶ�����һϵ�в����ͱ�����Ȼ��ִ����Դ��λ�ļ��㡣�����Ƕ���Ҫ���ֵ���ϸ���ͣ�

��ʼ���Ͳ�������:

theta_ �� varphi_�����������ռ�ĸ����Ǻ�ˮƽƫ��ǵķ�Χ���Զ�Ϊ��λ����Ȼ��ת��Ϊ���ȡ�
M �� N�����㸩���Ǻ�ˮƽƫ��ǵĳ��ȡ�
u, v, w�����������ռ���ÿ����ķ������ҡ�
sample_freq���������Ƶ�ʡ�
Audio_path��������Ƶ�ļ���·����
audioread����ȡ��Ƶ�ļ���
��˷�����:

mic_coordinate����������˷����е����ꡣ
��������:

search_coordinate������u, v, w���������ռ�����ꡣ
��������:

������chirp�źŵĲ�������chirp_length��chirp_interval��chirp_freq1 �� chirp_freq2��
�����γɲ����ͳ�ʼ��:

angle_cos��Ԥ�ȷ������ڴ洢�Ƕ����ҵľ���
����ʱ������ͼ:

ʹ��plot�������Ƶ�һ����Ƶͨ����ʱ�����С�
�����ӳٵ�:

����ÿ����˷������ÿ������������ӳٵ㡣
���ݴ���Ͳ����γ�:

ѡ��ο��ŵ��ķ�ֵ��Ϊͬ���㡣
����ÿ����ֵ����ȡ������˷�����ݶβ�ִ���ӳٺ���Ͳ����γɡ�
Energy_matrix���洢��ÿ�������ϵ�����ֵ��
������ӻ�:

ʹ��image��surf�������������������ά�����ֲ���
��Դ��λ:

�������������ҵ����ֵ��Ӧ�ĽǶȣ��Ӷ���λ��Դ��
�����ű���Ŀ�������ö����˷��¼����Ƶ��������λ��Դ����ͨ�����㲻ͬ�����ϵ����������ֲ����ҳ���Դ����ܵ�λ�á�
%}
theta_ = 0:180;
varphi_ = 0:180;
theta = theta_ * pi / 180;
varphi = varphi_ * pi /180;
M = length(theta);%������
N = length(varphi);%ˮƽƫ���
u = sin(theta)'*cos(varphi);
v = sin(theta)'*sin(varphi);
w = repmat(cos(theta)',1,N);
sample_freq = 44100;
Point = sample_freq * 0.6;
speed_sound = 343.0;
Audio_path = "C:\Users\82134\Desktop\";

%��˷������ȡ��wav�ļ�����ƺ����������ҵߵ���
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

data1 = T(:,1); % ��һ���ŵ�������
[ref_pk,ref_pos] = max(data1);
[pks, peak_loc_total] = findpeaks(data1,'MinPeakHeight',ref_pk*0.7,'MinPeakDistance',(chirp_interval -0.02)*fs);

% peak_loc_total = [peak_loc_total; ref_pos];

[n1, n2] = size(peak_loc_total); 
for beep_num = 1:n1
   waves = {};
    for data_num = 1:K
%         index_tmp = peak_loc_total(beep_num,1);%Ӧ��ѡ�ο��ŵ���peak��������ÿһ���ŵ���peak
        index_tmp = ref_pos; % ֻȥ������ǿ��һ���ź�
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
    set(gca,'ytick',[])  %��ȥy������ֵ
    set(gca,'xtick',[])  %��ȥx������ֵ
    set(gca,'Position',[ 0 0 1 1]);
    box off;
    grid on; 
    
    figure; %��άͼ
    surf(u,v,w,Energy_matrix, 'edgecolor', 'none');
end
