function combined_signal = delay_and_sum_beamforming_v5_yanzhi(waves,delay_point,start_beam,M,N,K)

rear = length(waves{1}(1,:)) - start_beam; 
sum_ = zeros(1,rear-start_beam+1);

E = zeros(M, N, length(sum_) ); 
%E= zeros(M,N);

%sum_ =  sum_ + waves{1}(1,start_beam:rear);
sum_ =  sum_ + waves{1}(1,start_beam : rear);

for j=1:M
    for k=1:N
        for i=2:K
            %sum_ = T(start_beam:rear,1) + T(start_beam + delay_point(j,k,i):rear + delay_point(j,k,i),i) + sum_;
            %sum_ = waves{1}(1,start_beam:rear) + waves{i}(1,start_beam + delay_point(j,k,i):rear + delay_point(j,k,i)) + sum_;
            sum_ = waves{i}(1,start_beam + delay_point(j,k,i):rear + delay_point(j,k,i)) + sum_;
        end
        sum_ = sum_ / 6; 
        %E(j,k) = sum(sum_.^2);
        E(j,k,:) = sum_; 
        
    end
end

combined_signal = E;


