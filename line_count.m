function [num,M_min,del_hang,del_lie]=line_count(M)
[a,~]=size(M);
num=0;
h=0;
del_hang=zeros(a,1);
del_lie=zeros(a,1);
for ii=1:a
    del=ii-h;
    [~,b]=size(find(M(del,:)==0));
    if   b>= 2
        M(del,:)=[];
        h=h+1;
        del_hang(ii)=ii;    %得到被覆盖的行数
        num=num+1;
    end
end
l=0;
for ii=1:a
    del=ii-l;
    [b,~]=size(find(M(:,del)==0));
    if  b >=1
        M(:,del)=[];
        l=l+1;
        del_lie(ii)=ii;    %得到被覆盖的列数
        num=num+1;
    end
end
M_min=min(min(M));
