function corrSeq = matched_filter_compute(data1, sampleSignal) 


%%%data initilization
sizeOfSample = length(sampleSignal); 
sizeOfData = length(data1); 
output = zeros(sizeOfData, 1); 


%%%Compute the correlation 
for i = 1  : sizeOfData-sizeOfSample+1
    output(i,1) = abs(dot(sampleSignal, data1(1, i:i+sizeOfSample-1))); 
    
end

%%%Return the final values
output = output/max(output); 
corrSeq = output; 
% 
% figure
% t1 = 1:length( output );
% t2 =  t1/44100;
% plot(t2, output,'LineWidth',2,'Color','b');
% hold on;


%%%Conduct the envelop detection 
window_size_for_envelope  = 10;  %%%The length of each frame


%%%Count the number of windows we need 
counter=0;
for i=1:(window_size_for_envelope-1):(length(output) - window_size_for_envelope)
    counter=counter+1;
end
envelope_m=zeros(counter,1);


%%%COmpute the values for each frame
counter=1;
for i=1:( window_size_for_envelope - 1 ):(length(output) - window_size_for_envelope)
    current_window = output(i:i+window_size_for_envelope-1,:);
    
    min_value=min(current_window);
    max_value=max(current_window);
    
    envelope_m(counter,1) = max_value; %%upper part
   % envelope_m(counter,1) = min_value; %%lower part
    
    counter=counter+1;
end


%%%Conduct the interpolation 
x=1:( window_size_for_envelope - 1 ):(length(output) - window_size_for_envelope);  %%%The position of the sampled points
y=envelope_m;  %%%The value of sampled points
x_i=1:(length(output) - window_size_for_envelope); %%%	

y_i = interp1(x,y,x_i,'spline');


%%%Return the results
corrSeq = [y_i  zeros( 1, (sizeOfData - length(y_i)) ) ];











