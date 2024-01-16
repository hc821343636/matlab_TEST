% 定义原始复矩阵的大小
n = 8; % 行数
m = 8; % 列数

% 生成一个随机复矩阵（示例复矩阵）

% 定义分解的矩阵数量
k = 8; % 分解成k个矩阵的乘积

% 初始化分解矩阵（实部和虚部）
A_real = rand(n, k);
A_imag = rand(n, k);
B_real = rand(k, k);
B_imag = rand(k, k);
C_real = rand(k, m);
C_imag = rand(k, m);
D_real = rand(k, m);
D_imag = rand(k, m);

% 设置迭代参数
max_iter = 1000; % 最大迭代次数
learning_rate = 0.2; % 学习率

% 迭代优化
for iter = 1:max_iter
    % 计算当前近似解
    approximated_matrix = (A_real + 1i * A_imag) * (B_real + 1i * B_imag) * (C_real + 1i * C_imag) * (D_real + 1i * D_imag);
    
    % 计算误差
    error = wnnk - approximated_matrix;
    
    % 计算梯度
    grad_A_real = -2 * real(error) * (B_real + 1i * B_imag)' * (C_real + 1i * C_imag)' * (D_real + 1i * D_imag)';
    grad_A_imag = -2 * imag(error) * (B_real + 1i * B_imag)' * (C_real + 1i * C_imag)' * (D_real + 1i * D_imag)';
    grad_B_real = -2 * (A_real + 1i * A_imag)' * real(error) * (C_real + 1i * C_imag)' * (D_real + 1i * D_imag)';
    grad_B_imag = -2 * (A_real + 1i * A_imag)' * imag(error) * (C_real + 1i * C_imag)' * (D_real + 1i * D_imag)';
    grad_C_real = -2 * (A_real + 1i * A_imag)' * (B_real + 1i * B_imag)' * real(error) * (D_real + 1i * D_imag)';
    grad_C_imag = -2 * (A_real + 1i * A_imag)' * (B_real + 1i * B_imag)' * imag(error) * (D_real + 1i * D_imag)';
    grad_D_real = -2 * (A_real + 1i * A_imag)' * (B_real + 1i * B_imag)' * (C_real + 1i * C_imag)' * real(error);
    grad_D_imag = -2 * (A_real + 1i * A_imag)' * (B_real + 1i * B_imag)' * (C_real + 1i * C_imag)' * imag(error);
    
    % 更新梯度时保留每行最大的两个数，其他设置为0
    [~, max_indices] = sort(A_real, 2, 'descend');
    for row = 1:n
        A_real(row, max_indices(row, 3:end)) = 0;
        A_imag(row, max_indices(row, 3:end)) = 0;
    end
    [~, max_indices] = sort(B_real, 1, 'descend');
    for col = 1:k
        B_real(max_indices(3:end, col), col) = 0;
        B_imag(max_indices(3:end, col), col) = 0;
    end
    [~, max_indices] = sort(C_real, 1, 'descend');
    for col = 1:k
        C_real(max_indices(3:end, col), col) = 0;
        C_imag(max_indices(3:end, col), col) = 0;
    end
    [~, max_indices] = sort(D_real, 1, 'descend');
    for col = 1:m
        D_real(max_indices(3:end, col), col) = 0;
        D_imag(max_indices(3:end, col), col) = 0;
    end
    
    % 计算RMSE
    rmse = norm(error, 'fro') / sqrt(n * m);
    
    fprintf('迭代 %d, RMSE: %f\n', iter, rmse);
    
    % 如果RMSE足够小，可以提前退出迭代
    if rmse < 1e-3
        break;
    end
end

% 显示结果
fprintf('原始复矩阵 original_matrix:\n');
disp(wnnk);

fprintf('分解的复矩阵 A (实部):\n');
disp(A_real);

fprintf('分解的复矩阵 A (虚部):\n');
disp(A_imag);

fprintf('分解的复矩阵 B (实部):\n');
disp(B_real);

fprintf('分解的复矩阵 B (虚部):\n');
disp(B_imag);

fprintf('分解的复矩阵 C (实部):\n');
disp(C_real);

fprintf('分解的复矩阵 C (虚部):\n');
disp(C_imag);

fprintf('分解的复矩阵 D (实部):\n');
disp(D_real);

fprintf('分解的复矩阵 D (虚部):\n');
disp(D_imag);

fprintf('近似分解的复矩阵 approximated_matrix:\n');
disp(approximated_matrix);
