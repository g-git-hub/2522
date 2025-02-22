function [best_params, best_le] = optimize_parameters(n_trials)
    if nargin < 1
        n_trials = 30;
    end
    
    % 初始化
    best_params = [10, 8/3, 31.49, 4.99];  % 基准参数
    best_le = [-Inf, -Inf, -Inf];
    best_score = -Inf;
    
    % 参数搜索范围（更保守）
    param_ranges = [
        9.8, 10.2;      % a: 10 ± 0.2
        2.6, 2.8;       % b: 8/3 ± 0.1
        31.2, 31.8;     % r: 31.49 ± 0.3
        4.8, 5.2        % c: 4.99 ± 0.2
    ];
    
    fprintf('开始参数优化 [%s]\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'));
    fprintf('搜索范围:\n');
    fprintf('a: [%.2f, %.2f]\n', param_ranges(1,:));
    fprintf('b: [%.2f, %.2f]\n', param_ranges(2,:));
    fprintf('r: [%.2f, %.2f]\n', param_ranges(3,:));
    fprintf('c: [%.2f, %.2f]\n', param_ranges(4,:));
    
    for i = 1:n_trials
        try
            % 生成参数（使用基准参数的小扰动）
            params = zeros(1, 4);
            for j = 1:4
                params(j) = param_ranges(j,1) + ...
                           (param_ranges(j,2) - param_ranges(j,1)) * rand();
            end
            
            % 计算Lyapunov指数
            [LE, ~] = compute_lyapunov_wolf(params);
            
            % 使用后半段数据计算平均值
            n_points = min(100, floor(size(LE,1)/2));
            if n_points < 10
                continue;
            end
            
            mean_le = mean(LE(end-n_points+1:end,:), 1);
            
            % 评估参数质量
            if ~any(isnan(mean_le)) && ~any(isinf(mean_le))
                % 混沌评分标准：
                % 1. 最大指数应为正
                % 2. 和应为负（耗散系统）
                % 3. 应有负指数（吸引子存在）
                score = mean_le(1) - 0.1*abs(sum(mean_le));
                
                if mean_le(1) > 0 && sum(mean_le) < 0 && ...
                   min(mean_le) < -0.1 && score > best_score
                    
                    best_score = score;
                    best_params = params;
                    best_le = mean_le;
                    
                    fprintf('\n发现更好的参数 (Trial %d):\n', i);
                    fprintf('参数: [%.4f, %.4f, %.4f, %.4f]\n', params);
                    fprintf('Lyapunov指数: [%.4f, %.4f, %.4f]\n', mean_le);
                    fprintf('得分: %.4f\n', score);
                end
            end
            
        catch ME
            fprintf('Trial %d: %s\n', i, ME.message);
            continue;
        end
    end
    
    if all(isinf(best_le))
        warning('未找到有效参数组合');
    end
end