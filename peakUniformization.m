%function [] = peakUniformization(peak_data)
max_Val=max(peak_data);
start=-1;
finish=-1;
for i=2:size(peak_data)-1
    if peak_data(i-1)==0 && peak_data(i)~=0
        start=i;
    end
    if peak_data(i)~=0&&peak_data(i+1)==0
        finish=i;
    end
    if start~=-1&&finish~=-1
        peak_data(start:finish)=peak_data(start:finish)+max_Val-max(peak_data(start:finish));
        start=-1;
        finish=-1;
    end
end
figure;
plot(x,peak_data);
ylim([max_Val-0.1,max_Val+0.15]);

    