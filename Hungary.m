function res=Hungary(N)
%输入的矩阵应N*N的
[a,~]=size(N);
%第一步每一行减去当前行最小值
for ii = 1:a
    N(ii,:)= N(ii,:)-min( N(ii,:));
end
%第二步每一列减去当前列最小值
for ii = 1:a
    N(:,ii)=  N(:,ii)-min( N(:,ii));
end
num=0;
while num~=a
    [num,N_min,del_hang,del_lie]=line_count(N);
    if num ~=a
        for ii=1:a
            if del_hang(ii)~=ii
                N(ii,:) =  N(ii,:)-N_min;
            end
            if del_lie(ii)==ii
            N(:,ii) =  N(:,ii)+N_min;
            end
        end
    else
        res=N;
    end
end
