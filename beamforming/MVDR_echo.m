clear; 
close all; 
clc;

theta_ = 0:180;
varphi_ = 0:180;
theta = theta_ * pi / 180;
varphi = varphi_ * pi /180;
M = length(theta);
N = length(varphi);
u = sin(theta)'*cos(varphi);
v = sin(theta)'*sin(varphi);
w = repmat(cos(theta)',1,N);
sample_freq = 44100;
Point = sample_freq*0.5;
speed_sound = 343.0;

Audio_path = "D:\matlab\matlabR2019b\bin\microphone_array\sound_test_4\office_44_1K_sampling_2-3K_Fre_5ms_duration_5s_Inter\1\";
f = 2530;

% [T(:,1),Fx] = audioread(Audio_path+"Audio Track.wav",[offset Point+offset]);
% [T(:,2),Fx] = audioread(Audio_path+"Audio Track-2.wav",[offset Point+offset]);
% [T(:,3),Fx] = audioread(Audio_path+"Audio Track-3.wav",[offset Point+offset]);
% [T(:,4),Fx] = audioread(Audio_path+"Audio Track-4.wav",[offset Point+offset]);
% [T(:,5),Fx] = audioread(Audio_path+"Audio Track-5.wav",[offset Point+offset]);
% [T(:,6),Fx] = audioread(Audio_path+"Audio Track-6.wav",[offset Point+offset]);
[T(:,1),fs] = audioread(Audio_path+"Audio Track-4.wav");
[T(:,2),fs] = audioread(Audio_path+"Audio Track-5.wav");
[T(:,3),fs] = audioread(Audio_path+"Audio Track-6.wav");
[T(:,4),fs] = audioread(Audio_path+"Audio Track.wav");
[T(:,5),fs] = audioread(Audio_path+"Audio Track-2.wav");
[T(:,6),fs] = audioread(Audio_path+"Audio Track-3.wav");

%%%Parameters for chirp signals
chirp_length = 0.002;  %%%the length of chirp 2 ms
chirp_interval = 0.5;  %%%the interval of chirp
chirp_freq1 = 2000; 
chirp_freq2 = 3000; 
win_beam = round((0.04)*fs);
start_beam = round(0.0035*fs);


mic_coordinate(1,:) = [0      0   0]; %1# microphone
mic_coordinate(2,:) = [0.05   0   0];
mic_coordinate(3,:) = [0.075  0   0.0425];
mic_coordinate(4,:) = [0.05   0   0.085];
mic_coordinate(5,:) = [0      0   0.085];
mic_coordinate(6,:) = [-0.025 0   0.0425];   
% mic_coordinate = -mic_coordinate;
search_coordinate(:,:,1)=u;
search_coordinate(:,:,2)=v;
search_coordinate(:,:,3)=w;
K = size(mic_coordinate,1);
angle_cos = zeros(M,N,K);
% Calculate delay points
for i=2:K
    for j=1:M
        for k=1:N
            angle_cos(j,k,i)=squeeze(mic_coordinate(i,:))*squeeze(search_coordinate(j,k,:));
        end
    end
end

peak_loc_total = [];
for i = 1:K
   data1 = T(:,i);
   [pks, peak_loc] = findpeaks(data1,'MinPeakHeight',0.05,'MinPeakDistance',(chirp_interval -0.02)*fs);
   peak_loc_total = [peak_loc_total peak_loc]; 
end

[n1, n2] = size(peak_loc_total); 
count = 0;
for chirp_num = 5:5
    waves = {};
    index_tmp = peak_loc_total(chirp_num, 1);
    for data_num = 1:n2
       data_seg = T(index_tmp: index_tmp+win_beam,data_num);
       TT(:,data_num) = data_seg'; 
    end
    Y = fft(TT);
%     P1 = abs(Y/Point);
%     P2 = P1(1:Point/2+1,:);
%     P2(2:end-1,:) = 2*P2(2:end-1,:);

    S2 = Y(f*Point/sample_freq+1,:).';
    R = S2*S2'*10^(30/10)+eye(6);
    E = zeros(M,N);
    M_ = [0:5];
    for j=1:M
        for k=1:N
            a = exp(1i*2*pi*squeeze(angle_cos(j,k,:))*f/speed_sound);
            W = R\a / (a'/R*a);
            E(j,k) = W'*R*W;
        end
    end
end

    [x,y]=find(E==max(max(E)));
    theta_max = theta_(x);
    varphi_max = varphi_(y);
    fprintf("theta: %d, varphi: %d \n",theta_max,varphi_max);
    fprintf("*****************************\n");

E = abs(E);
figure;
surf(u,v,w,E, 'edgecolor', 'none');
axis('square')

figure;
image(E,'CDataMapping','scaled');
colorbar('off');
%caxis([0,0.2]);
set(gca,'ytick',[])  %隐去y轴坐标值
set(gca,'xtick',[])  %隐去x轴坐标值
box off;
grid on;
