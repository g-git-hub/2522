function [best_results, all_results] = GridSearch()
    templates = SystemTemplates();
    num_systems = length(templates);
    
    % 参数网格定义
    param_grids = {
        % 系统1的参数网格
        {
            linspace(20, 40, 10),  % a
            linspace(2, 5, 10),    % b
            linspace(20, 45, 10),  % r
            linspace(0.5, 2, 10),  % c
            linspace(0.5, 2, 10),  % d
            linspace(0.5, 2, 10)   % k
        },
        % 系统2的参数网格
        {
            linspace(8, 15, 10),   % sigma
            linspace(2, 4, 10),    % beta
            linspace(25, 40, 10),  % rho
            linspace(0.5, 2, 10),  % c
            linspace(0.5, 2, 10),  % d
            linspace(0.5, 2, 10)   % k
        },
        % ... 其他系统的参数网格
    };
    
    best_results = cell(num_systems, 1);
    all_results = cell(num_systems, 1);
    
    parfor sys_idx = 1:num_systems
        system = templates{sys_idx};
        grid = param_grids{sys_idx};
        
        fprintf('开始搜索系统 %d...\n', sys_idx);
        
        results = struct('params', [], 'le', [], 'score', []);
        count = 0;
        
        % 网格搜索
        for p1 = grid{1}
            for p2 = grid{2}
                for p3 = grid{3}
                    for p4 = grid{4}
                        for p5 = grid{5}
                            for p6 = grid{6}
                                params = [p1,p2,p3,p4,p5,p6];
                                try
                                    [LE, ~] = compute_lyapunov_wolf(@(t,X) system(t,X,params));
                                    mean_le = mean(LE(end-20:end,:), 1);
                                    
                                    if ~any(isnan(mean_le)) && ~any(isinf(mean_le))
                                        count = count + 1;
                                        results(count).params = params;
                                        results(count).le = mean_le;
                                        
                                        % 评分函数
                                        score = EvaluateLE(mean_le);
                                        results(count).score = score;
                                        
                                        if mod(count, 100) == 0
                                            fprintf('系统%d: 已测试%d组参数\n', sys_idx, count);
                                        end
                                    end
                                catch
                                    continue;
                                end
                            end
                        end
                    end
                end
            end
        end
        
        % 对结果进行排序
        [~, idx] = sort([results.score], 'descend');
        best_results{sys_idx} = results(idx(1));
        all_results{sys_idx} = results(idx(1:min(100,length(idx))));  % 保存前100个结果
    end
end

function score = EvaluateLE(le)
    % 改进的评分函数
    if le(1) <= 0
        score = -inf;
        return;
    end
    
    score = le(1) * 5 + ...                    % 最大LE权重
           (le(1) > 0.5) * 3 + ...            % 奖励较大LE
           (le(1) > 1.0) * 5 + ...            % 奖励很大LE
           (le(1) > 1.5) * 8 + ...            % 奖励极大LE
           (le(2) > 0) * le(2) * 2 - ...      % 奖励超混沌
           0.1 * abs(sum(le));                 % 保持适度耗散
           
    % 稳定性检查
    if any(abs(le) > 10)
        score = -inf;
    end
end