 %%  ��������  %%%%%%%%%%%
fs = 44100;             % ����Ƶ��
bps = 16;               % λ�� bit
Tfile =10;              % ����Ƭ�ε���ʱ��,��λs
vol = 0;                % ����Ƭ�ε����� 0db
filename = ('2ms_44100_2-3KHz_60s_chirp.wav'); %����ļ����ƣ��Լ����壬��׺��.wav��
%% ����ʱ������
Df = 1;                 % Ƶ�ʼ����Ĭ��1
T = 1/fs;               % ��������
N = fs/Df;              % ���е���
time = Tfile*(N-1)*T;   % �ļ���ʱ��
repeat = 120;
%t = 0:T:(0.5+duration)*repeat;           % ����ÿһ���������Ӧ��ʱ������
%% ������������
duration = 0.002; %ÿ��chirp����ʱ��
feq1 = 2000;
feq2 = 3000;
tt = 0:1/fs:duration;

x = chirp(tt,feq1,duration,feq2,'linear');
h = hamming(length(x));
x = x.*h';
figure;
plot(tt(1:length(x)),x,'LineWidth',3,'Color','b');
set(gca,'ytick',[])  %��ȥy������ֵ
set(gca,'xtick',[])  %��ȥx������ֵ

y=[x,zeros(1,0.5*fs)];%ÿ�μ��0.5��
figure;
spectrogram(y,256,250,256,fs,'yaxis');
t = 0:T:(length(y)/fs);
figure;
plot(t(1:length(y)),y);


output = repmat(y,1,repeat); % chirp�ظ�����
t = 0:T:(length(output)/fs);
%y = chirp(mod(t,0.005),11000,T,14000);
%fb = (1+square(4*pi*t,1))*0.5;
%yn = y .* fb;
%f0 = 20;
%fe = 200;
%ft=f0+(fe-f0)*mod(t,0.02);
%y0 = abs(fft(y));
%f = 0:fs/length(y0);
% figure;      
% spectrogram(output,256,250,256,fs,'yaxis');

output =  bandpass(output,[1500 3500],fs);
figure;      
spectrogram(output,256,250,256,fs,'yaxis');

figure;
plot(t(1:length(output)),output)
hold on;

%% д���ļ�
audiowrite(filename,output,fs,'BitsPerSample',bps);  % �洢.wav��Ƶ�ļ����ļ����ڲ�������������
%audiowrite('2ms_chirp_soundv2_2-3k_44100.wav',x,fs,'BitsPerSample',bps);