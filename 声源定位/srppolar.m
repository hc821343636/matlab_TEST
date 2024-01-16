function [finalpos,finalsrp,finalfe]=srppolar(s, mic_loc, fs, lsb, usb)
%%
%{
函数定义与初始化:
function [finalpos,finalsrp,finalfe]=srppolar(s, mic_loc, fs, lsb, usb)：定义了一个名为srppolar的函数，其输入参数分别为信号s，麦克风位置mic_loc，采样频率fs，以及两个可选参数lsb和usb。返回值为finalpos, finalsrp, finalfe。
参数检查和默认值设置:

if nargin < 5, usb=[2 0 6]; end等：如果用户没有提供所有的参数，这些行会设置默认值。
初始化和准备工作:

L = size(s,1);：获取信号帧的长度。
M = size(mic_loc,1);：计算麦克风的数量。
np=M*(M-1)/2;：计算独立麦克风对的数量。
声速和相关常数计算:

计算声速和一些与采样频率和声速相关的常数，如magiconst。
麦克风之间的最大距离计算:

使用pdist函数计算麦克风位置之间的距离，并找出最大值。
广义互相关(GCC-PHAT)处理:

对信号进行FFT变换。
计算所有独立麦克风对的互相关，并应用PHAT权重。
将得到的互相关结果存储在yv中。
对GCC-PHAT进行立方样条插值:

提高GCC-PHAT的分辨率。
声源定位算法的核心部分:

初始化并进行声源定位搜索。
通过计算不同方向上的SRP-PHAT值，确定声源的大致方向。
绘制极坐标图:

使用MATLAB的绘图功能展示声源定位的结果。
函数fe的定义：

fe函数用于评估搜索空间中点的SRP-PHAT值，并返回该点的位置和对应的SRP-PHAT值。
%}
warning off all


if nargin < 5, usb=[2 0 6]; end
if nargin < 4, lsb=[-2 -1 0]; end
if nargin < 3, fs=20000; end


L = size(s,1); %%% determine frame-length
M = size(mic_loc,1); %%% number of microphones (also number of channels of the input X).
np=M*(M-1)/2; %%%number of independent pairs

%steplength = L/4;  %%% non-overlapped length of a frame (75% overlap)
dftsize = L;       %%% dft size equals to frame size
temperatureC=24.0;
speedofsound=331.4*sqrt(1.0+(temperatureC/273));
magiconst=10*fs/speedofsound;  

%% Determine the maximum end-fire length (in samples) of a microphone pairs:
mdist=pdist(mic_loc); %pdist获取两点间距离
efmax=max(mdist);
efs=2*(fix(efmax*fs/speedofsound)); %%% End-fire is symmetric about the 0th sample so multiplying by 2.
%efs=801;
hefs=round(efs/2); %%%half of the end-fire region
%% Get the linear-indices of 'np' independent mic-pairs:
% w=1:M;
% wn=[0:M-1]'*M;
% fm1=repmat(wn,1,M);
% fm2=repmat(w,M,1);
% mm=fm1+fm2;
% tr=triu(mm,1);%%Get the upper half above the main diagonal.
% gidM=nonzeros(tr'); %%%keep only non-zero values -> linear indices of the pairs
% clear fm1 fm2 w wn

%% Doing the GCC-PHAT:


sf=fft(s,dftsize);                    %%%FFT of the original signals  
csf=conj(sf);
yv=zeros(np,efs);                     %%%%Initialize yv to store the SRP-PHAT samples

p=1;
for i=1:M-1
      su1mic=sf(:,i)*ones(1,M);      %%%Create M copies of each signal's FFT
      prodall=su1mic.*csf;           %%%%Calculate the cross-power spectrum: = fft(x1).*conj(fft(x2))
      ss=prodall(:,i+1:M);           %%%% ss will be the cross-power spectra of microphone pairs (i,i+1), (i,i+2)...(i,M)
      ss=ss./(abs(ss)+eps);          %%%% PHAT weighting
      
      ssifft=real(ifft(ss,dftsize)); %%%% Get the GCC-PHAT ssifft, which is the IFFT of the cross-power spectra
      %newssifft=[ssifft(end-hefs:end,:); ssifft(1:efs-hefs-1,:)]; %%% Only select 'efs' samples (the beginning+end portions)
      newssifft=[ssifft(end-hefs+1:end,:); ssifft(1:efs-hefs,:)];
      newssifftr = newssifft';          %%%% Transpose it
      yv(p:p+M-i-1,:)=newssifftr;    %%%%Store in yv
      p=p+M-i;                       %%%% Update the current index of yv
      clear su1mic prodall ss ssifft newssifft newssifftr
end

%% Doing cubic-splines interpolation (factor of 10) on the GCC-PHAT:

xx=[1:.1:efs];
x=[1:efs];
yintp=spline(x,yv,xx);
yintpt=yintp';

%% Initialize to do SRC:

efsintp=length(xx)/2;

row1=([0:np-1]*2*efsintp)'; 

randpts=3000;  %%% J0 in SRC. Depending on the size of the search volume, choose an appropriate value here (Here, 3000 is for a V_{search}= 4m x 1m x 6m) 
npoints=100;   %%% Best N points. Again, choose an appropriate number according to your problem.

varphi = linspace(0, pi, 181);
v_length = length(varphi);
search_r = 1.0;

srp_value = [];
for i= 1:v_length
   src_pos = [search_r*cos(varphi(i)) search_r*sin(varphi(i)) 0.5] ;
   [yval1,~] = fe(src_pos,magiconst,mic_loc,yintpt,row1,efsintp);
   srp_value = [srp_value yval1];
end

% 创建一个 figure
figure;
ax = polaraxes;
polarplot(varphi, srp_value, '-r', 'LineWidth', 2);
title('Polar Plot Example');

% 可选：如果需要将0度方向放在极坐标图的上方
ax = gca;
ax.ThetaZeroLocation = 'top';

% 显示刻度标签
rticks([-1 -0.5 0 0.5 1]); % 修改刻度位置
rticklabels({'-1', '-0.5', '0', '0.5', '1'}); % 修改刻度标签

% 可选：添加网格线
grid on;

finalpos=0;
finalsrp=0;
finalfe=0; 

%% Functional Evaluation sub-routine (calculate SRP-PHAT value for a point in the search space):

    function [yval1,position1] = fe(src_pos,magiconst,mic_loc,yintpt,row1,efsintp)
        %%%This function evaluates an 'fe' in the search space, gives back the 'fe' position and its SRP-PHAT value.


        M=size(mic_loc,1); %%%number of mics

        %%%generate a random point in the search space:
        x=src_pos;   %This vector x defines the coordinates (x,y,z) of the rand point x

        %%%Find the distances from M-microphones to the point:
        a1=ones(M,1);
        xx1=a1*x;
        xdiff1=xx1-mic_loc;
        dists=sqrt(sum(xdiff1.*xdiff1,2));


        %%%%Differences in distances:
        ddiffs_ones=ones(M,1)*dists';
        ddm=ddiffs_ones-ddiffs_ones';

        %%% Calculate the TDOA index:
      %  v=ddm(gidM);
        v=nonzeros(tril(ddm,0));

        %ddiffsi32=int32(round(magiconst*v+efsintp));
        %ddiffsi32=floor(magiconst*v+efsintp)+1; 
        ddiffsi32=round(magiconst*v+efsintp+1.5); 
        row=row1+ddiffsi32;   
        
        %%%Pull out the GCC-PHAT value corresponding to that TDOA:
        v1=yintpt(row);
        %%% SRP-PHAT value of the point:
        yval1=sum(v1);
        position1=x; 
    end
end