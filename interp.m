% 生成示例信号
x = linspace(0, 10, 101);
y = sin(x);
y(x>5 & x<7) = 0;
plot(x,y,'.-')

% 对连续0信号进行插值
idx = find(y == 0);
for ii = 1:length(idx)
    if ii == 1 || idx(ii) ~= idx(ii-1) + 1 % 判断是否是连续0的序列的开头
        start_idx = idx(ii);
        end_idx = start_idx;
        while end_idx < length(y) && y(end_idx+1) == 0 % 找到连续0序列的末尾
            end_idx = end_idx + 1;
        end
        if end_idx == length(y) % 特判连续0序列在末尾的情况
            end_idx = end_idx + 1;
        end
        x_interp = x(start_idx:end_idx); % 进行插值的样本点
        y_interp = y(start_idx:end_idx);
        y_interp(y_interp == 0) = nan; % 将0替换为nan，以便于interp1函数的使用
        y_interp = interp1(x_interp(~isnan(y_interp)), y_interp(~isnan(y_interp)), x_interp, 'linear');
        y(start_idx:end_idx) = y_interp;
    end
end

hold on
plot(x,y,'.-')
legend('原始信号','插值后信号')
