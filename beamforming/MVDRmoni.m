%单目标MVDR和CBF波束形成，通过调整SNR改变信噪比，根据实验所需信号来选择性注释
clear;clc;fclose('all');
%声源、水平阵、环境参数设置
%环境参数
c = 343;   %声速1500m/s
SNR = 20;   %信噪比为20dB
%声源（假定在很远处的点源，近场球面波,到达很远处的水平阵的时候可近似为平面波）
SL = 140;   %信号能量140dB
r = 7000;   %与水平阵首阵元相距7000m
f = 50;     %信号频率50Hz
lambda=c/f; %波长
angle=30;   %入射角
%水平阵
M = 9;             %水平阵阵元数9
d = c/(2*f);       %阵元间距15m
angles=-90:0.1:90; %检测角度范围
Nsnapshot=15;      %快拍数（结合多次发射信号）
%计算各阵元接收信号x=s+n,即信号加噪声
%阵列响应向量
v=sqrt(M)\exp(-1j*2*pi*(d*sin(angle*pi/180)/lambda)*(-(M-1)/2:(M-1)/2)'); 
%发射信号
s=sqrt(10^(SL/10))*exp(1j*2*pi*f*(1:Nsnapshot));
%高斯白噪声
n=sqrt(10^((SL-SNR)/10))*(randn(M,Nsnapshot)+1j*randn(M,Nsnapshot))/sqrt(2);
%接收信号
x=sqrt(M)*v*s+n;
%计算响应向量和波束形成响应
%Rx=(M*v*(s*s')*v'+10^((SL-SNR)/10)*eye(M))/Nsnapshot;%理想信号
Rx=(x*x')/Nsnapshot;%加入高斯白噪声的真实信号
%驾驶向量
c=sqrt(M)\exp(-1j*2*pi*(-(M-1)/2:(M-1)/2)'*(d*sin(angles*pi/180)/lambda));
%直接求解闭式解
Cmv1=(Rx\c)/(diag(diag((c'/Rx)*c)));
y1=diag(Cmv1'*Rx*Cmv1);     %MVDR响应
%y1=diag(c'*Rx*c);     %CBF响应
y1=abs(y1)/max(abs(y1));
%绘图
figure(1);
plot(angles,10*log10(y1),'k','linewidth',2);
xlabel('Angle(deg)','Fontsize',15);ylabel('Power Response(dB)','Fontsize',15) ;
title('单目标MVDR波束形成（实际信号）','Fontsize',20);
