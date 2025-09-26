% 设定参数
N = 10; M = 5; k = 3;
Q = (randn(M,N) + 1i*randn(M,N)) / sqrt(2);

% 随机生成 k-sparse 信号
b_true = zeros(N,1);
support_true = randperm(N,k);
b_true(support_true) = randn(k,1) + 1i*randn(k,1);

% 生成观测
A = Q * b_true;

supportest = bruteforce(A,Q,k);

% % 初始化
% support_list = nchoosek(1:N, k);
% num_supports = size(support_list, 1);
% min_residual = inf;
% 
% % 穷举所有支持集
% for i = 1:num_supports
%     S = support_list(i,:);
%     Qs = Q(:, S);
%     b_s = Qs \ A;
%     residual = norm(A - Qs * b_s);
% 
%     if residual < min_residual
%         min_residual = residual;
%         best_support = S;
%         best_b_s = b_s;
%     end
% end
% 
% % 构造完整的估计向量
% b_est = zeros(N,1);
% b_est(best_support) = best_b_s;
% 
% % 输出结果
% disp('真实支持集:');
% disp(sort(support_true));
% disp('估计支持集:');
% disp(sort(best_support));
% disp(['重建误差: ', num2str(norm(b_est - b_true))]);
