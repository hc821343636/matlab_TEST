function [place,res]=linear_assignment(M,N)
%N是n维矩阵,N是经过Hungary处理的
%M是未处理前的
[a,~]=size(N);
x=0;
place=zeros(1,a);
res=zeros(1,a);
judge=zeros(1,a);
while find(N==0)
    for ii=1:a
    judge(ii)=length(find(N(ii,:)==0));
    end
    judge(find(judge==0))=[];
    if min(judge)==1
     for ii=1:a
        if length(find(N(ii,:)==0))==1     %先选出行中只有1个0
            x=x+1;
            place(x)=ii+(find(N(ii,:)==0)-1)*a; %得到矩阵中的位置
            h=find(N(ii,:)==0);
            N(ii,:)=1./zeros(1,a);
            N(:,h)=1./zeros(a,1);
        end
     end
    end
    
    for ii=1:a
    judge(ii)=length(find(N(ii,:)==0));
    end
    judge(find(judge==0))=[];
    
    if min(judge)==2
       x=x+1;
    q=find(N==0);
    place(x)=q(1);
    N(mod(q(1),a),:)=1./zeros(1,a);
    N(:,fix(q(1)/a)+1)=1./zeros(a,1);  
    end
end
[place,~]=sort(place);
for ii=1:length(place)
    res(ii)=M(place(ii));
end
 
 
