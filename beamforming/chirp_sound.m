 %%  参数定义  %%%%%%%%%%%
fs = 44100;             % 采样频率
bps = 16;               % 位深 bit
Tfile =10;              % 声音片段的总时长,单位s
vol = 0;                % 声音片段的音量 0db
filename = ('2ms_44100_2-3KHz_60s_chirp.wav'); %输出文件名称（自己定义，后缀是.wav）
%% 生成时间序列
Df = 1;                 % 频率间隔，默认1
T = 1/fs;               % 采样周期
N = fs/Df;              % 序列点数
time = Tfile*(N-1)*T;   % 文件总时长
repeat = 120;
%t = 0:T:(0.5+duration)*repeat;           % 生成每一个采样点对应的时间序列
%% 生成声道波形
duration = 0.002; %每次chirp持续时间
feq1 = 2000;
feq2 = 3000;
tt = 0:1/fs:duration;

x = chirp(tt,feq1,duration,feq2,'linear');
h = hamming(length(x));
x = x.*h';
figure;
plot(tt(1:length(x)),x,'LineWidth',3,'Color','b');
set(gca,'ytick',[])  %隐去y轴坐标值
set(gca,'xtick',[])  %隐去x轴坐标值

y=[x,zeros(1,0.5*fs)];%每次间隔0.5秒
figure;
spectrogram(y,256,250,256,fs,'yaxis');
t = 0:T:(length(y)/fs);
figure;
plot(t(1:length(y)),y);


output = repmat(y,1,repeat); % chirp重复次数
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

%% 写入文件
audiowrite(filename,output,fs,'BitsPerSample',bps);  % 存储.wav音频文件，文件名在参数定义里设置
%audiowrite('2ms_chirp_soundv2_2-3k_44100.wav',x,fs,'BitsPerSample',bps);