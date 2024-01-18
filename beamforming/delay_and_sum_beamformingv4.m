function Energy_matrix = delay_and_sum_beamformingv4(waves,delay_point,start_beam,M,N,K)
rear = length(waves{1}(1,:)) - start_beam;
E= zeros(M,N);
sum_ = zeros(1,rear-start_beam+1);
sum_ =  sum_ + waves{1}(1,start_beam:rear);
for j=1:M
    for k=1:N
        for i=2:K
            %sum_ = T(start_beam:rear,1) + T(start_beam + delay_point(j,k,i):rear + delay_point(j,k,i),i) + sum_;
            %sum_ = waves{1}(1,start_beam:rear) + waves{i}(1,start_beam + delay_point(j,k,i):rear + delay_point(j,k,i)) + sum_;
            sum_ = waves{i}(1,start_beam + delay_point(j,k,i):rear + delay_point(j,k,i)) + sum_;
        end
        
        sum_ = sum_ / 6;
%         if j<=20 ||j>=160 ||k<=20 ||k>=160
%             sum_ = sum_ * 0.8;
%         end
        E(j,k) = sum(sum_.^2);
    end
end

Energy_matrix = E;

end

