function [r] = target_fun(B,Fn,A,P,D)
len=size(A,1);
F_hat=eye(len,len);
F_hat=P*F_hat;
F_hat=F_hat*A(:,:,1);
F_hat=F_hat*D;
F_hat=F_hat*A(:,:,3)*A(:,:,2)*A(:,:,1);
%r=sqrt(mean((B*Fn-F_hat).^2));


errorMatrix = F_hat - B*Fn;
% 计算误差平方的均值
r=sqrt(sum(sum(errorMatrix .* errorMatrix)))/len;
disp(sum(errorMatrix .* errorMatrix));

disp(norm(errorMatrix, 'fro') / sqrt(len *len));

end

