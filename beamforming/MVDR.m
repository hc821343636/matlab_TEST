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
offset = 1000;

f = 600;

[T(:,1),Fx] = audioread(Audio_path+"Audio Track.wav",[offset Point+offset]);
[T(:,2),Fx] = audioread(Audio_path+"Audio Track-2.wav",[offset Point+offset]);
[T(:,3),Fx] = audioread(Audio_path+"Audio Track-3.wav",[offset Point+offset]);
[T(:,4),Fx] = audioread(Audio_path+"Audio Track-4.wav",[offset Point+offset]);
[T(:,5),Fx] = audioread(Audio_path+"Audio Track-5.wav",[offset Point+offset]);
[T(:,6),Fx] = audioread(Audio_path+"Audio Track-6.wav",[offset Point+offset]);

Y = fft(T);
P1 = abs(Y/Point);
P2 = P1(1:Point/2+1,:);
P2(2:end-1,:) = 2*P2(2:end-1,:);

S2 = Y(f*Point/sample_freq+1,:).';
R = S2*S2'*10^(30/10)+eye(6);
% Pyy = [1 : 6];
% for i = 1 : 6
%     Pyy(i) =phase(S2(i));               %ËÆ°ÁÆóÁõ∏‰Ωç
%     Pyy(i) = Pyy(i) * 180 /pi;          %Êç¢ÁÆó‰∏∫ËßíÂ∫?
% end
% plot([1:6],Pyy);

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

E = zeros(M,N);
M_ = [0:5];
for j=1:M
    for k=1:N
        a = exp(1i*2*pi*squeeze(angle_cos(j,k,:))*f/speed_sound);
        W = R\a / (a'/R*a);
        E(j,k) = W'*R*W;
    end
end

E = abs(E);
figure;
surf(u,v,w,E, 'edgecolor', 'none');
axis('square')

figure;
image(E,'CDataMapping','scaled');
colorbar('off');
%caxis([0,0.2]);
set(gca,'ytick',[])  %“˛»•y÷·◊¯±Í÷µ
set(gca,'xtick',[])  %“˛»•x÷·◊¯±Í÷µ
box off;
grid on;
