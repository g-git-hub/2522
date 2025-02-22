% 主程序：增强型混沌系统分析
clear all;
close all;

% 1. 寻找最优参数
[best_params, best_le] = optimize_parameters();

% 2. 显示优化结果
fprintf('\n最终系统参数:\n');
param_names = {'a (耦合强度)', 'b (阻尼系数)', 'r (主控参数)', ...
               'c (非线性增益)', 'd (二次项系数)', 'k (三角函数调制)', ...
               'm (双曲函数调制)', 'n (指数项系数)'};
for i = 1:length(best_params)
    fprintf('%s = %.4f\n', param_names{i}, best_params(i));
end

fprintf('\nLyapunov指数谱:\n');
fprintf('λ1 = %.4f\n', best_le(1));
fprintf('λ2 = %.4f\n', best_le(2));
fprintf('λ3 = %.4f\n', best_le(3));

% 3. 计算Kaplan-Yorke维数
if sum(best_le(1:2)) > 0 && best_le(3) < 0
    ky_dim = 2 + sum(best_le(1:2))/abs(best_le(3));
    fprintf('\nKaplan-Yorke维数: %.4f\n', ky_dim);
end

% 4. 进行可视化分析
visualization_analysis(best_params);