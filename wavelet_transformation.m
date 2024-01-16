function [xrec] = wavelet_transformation(data)

data_row = data';
[c,l] = wavedec(data_row,3,'db5');
thr = wthrmngr('dw1dcompGBL','bal_sn',c,l);
c(l(1)+1:end) = c(l(1)+1:end).*(c(l(1)+1:end)>thr);
xrec = waverec(c,l,'db5');


%{
figure;
L = 5;
data_row=[data_row,data_row(end-3:end)];
swc = swt(data_row,L,'haar');
swcnew = swc;
ThreshML = wthrmngr('sw1ddenoLVL','penallo',swc,3);
for jj = 1:L
    swcnew(jj,:) = wthresh(swc(jj,:),'h',ThreshML(jj));
end
noisbloc_denoised = iswt(swcnew,'haar');
plot(data_row)
hold on
plot(noisbloc_denoised,'r','linewidth',1)
legend('Original','Denoised')
%}