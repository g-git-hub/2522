function [best_params, best_le] = optimize_parameters()
    % 基本参数设置
    num_particles = 50;
    max_iterations = 50;
    max_time_per_iter = 60; % 每次迭代最大允许时间(秒)
    
    % 参数范围
    param_ranges = [
        [9, 12];    % a
        [2.5, 3.0]; % b
        [30, 35];   % r
        [3.5, 5.0]  % c
    ];
    
    % 初始化
    num_params = size(param_ranges, 1);
    particles = zeros(num_particles, num_params);
    velocities = zeros(num_particles, num_params);
    p_best = zeros(num_particles, num_params);
    p_best_scores = -inf(num_particles, 1);
    g_best = zeros(1, num_params);
    g_best_score = -inf;
    
    % 随机初始化粒子位置
    for i = 1:num_particles
        particles(i,:) = param_ranges(:,1)' + ...
            rand(1,num_params) .* (param_ranges(:,2)' - param_ranges(:,1)');
    end
    
    % 初始化速度
    velocities = randn(num_particles, num_params) * 0.05;
    
    % PSO参数
    w = 0.7;  % 惯性权重
    c1 = 1.5; % 个体学习因子
    c2 = 1.5; % 群体学习因子
    
    fprintf('开始优化...\n');
    
    % 主循环
    for iter = 1:max_iterations
        tic; % 开始计时
        fprintf('迭代 %d: 开始计算\n', iter);
        
        % 创建临时数组存储结果
        temp_scores = zeros(num_particles, 1);
        temp_LE = zeros(num_particles, 3);
        valid_results = false(num_particles, 1);
        
        % 并行计算每个粒子
        parfor i = 1:num_particles
            try
                [LE, ~] = compute_lyapunov_wolf(particles(i,:));
                mean_le = mean(LE(end-10:end,:), 1);
                
                if ~any(isnan(mean_le)) && ~any(isinf(mean_le))
                    temp_LE(i,:) = mean_le;
                    temp_scores(i) = calculate_score(mean_le);
                    valid_results(i) = true;
                end
            catch
                continue;
            end
        end
        
        % 更新个体最优
        for i = 1:num_particles
            if valid_results(i) && temp_scores(i) > p_best_scores(i)
                p_best_scores(i) = temp_scores(i);
                p_best(i,:) = particles(i,:);
            end
        end
        
        % 更新全局最优
        [max_score, max_idx] = max(p_best_scores);
        if max_score > g_best_score && valid_results(max_idx)
            g_best_score = max_score;
            g_best = p_best(max_idx,:);
            
            % 输出当前最优结果
            ky_dim = 2 + (temp_LE(max_idx,1) + temp_LE(max_idx,2))/abs(temp_LE(max_idx,3));
            fprintf('\n迭代 %d:\n', iter);
            fprintf('参数值:\n');
            fprintf('a = %.4f\n', g_best(1));
            fprintf('b = %.4f\n', g_best(2));
            fprintf('r = %.4f\n', g_best(3));
            fprintf('c = %.4f\n', g_best(4));
            fprintf('Lyapunov指数: [%.4f %.4f %.4f]\n', temp_LE(max_idx,:));
            fprintf('Kaplan-Yorke维数: %.4f\n', ky_dim);
        end
        
        % 更新粒子位置和速度
        r1 = rand(num_particles, num_params);
        r2 = rand(num_particles, num_params);
        
        velocities = w * velocities + ...
                    c1 * r1 .* (p_best - particles) + ...
                    c2 * r2 .* (repmat(g_best, num_particles, 1) - particles);
                
        % 限制速度
        max_vel = 0.1 * (param_ranges(:,2) - param_ranges(:,1))';
        velocities = min(max(velocities, -max_vel), max_vel);
        
        % 更新位置
        particles = particles + velocities;
        
        % 确保在参数范围内
        for i = 1:num_particles
            particles(i,:) = min(max(particles(i,:), param_ranges(:,1)'), param_ranges(:,2)');
        end
        
        % 检查迭代时间
        iter_time = toc;
        fprintf('迭代 %d: 完成，用时 %.2f 秒\n', iter, iter_time);
        
        % 如果单次迭代超时，提前结束
        if iter_time > max_time_per_iter
            fprintf('迭代时间超过限制，优化提前结束\n');
            break;
        end
    end
    
    % 返回最优结果
    best_params = g_best;
    [LE, ~] = compute_lyapunov_wolf(best_params);
    best_le = mean(LE(end-10:end,:), 1);
end

function score = calculate_score(le)
    % 简单的评分函数，主要关注混沌特性
    if le(1) > 0 && le(2) < 0.1 && le(2) > -0.1
        score = le(1);  % 正的最大Lyapunov指数
    else
        score = -inf;
    end
end