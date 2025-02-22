function [LE, t_span] = compute_lyapunov(params, tmax, dt)
    if nargin < 2
        tmax = 20;    % 减少总时间避免不稳定性
    end
    if nargin < 3
        dt = 0.01;    % 增加时间步长
    end
    
    % 系统维数
    n = 3;
    
    % 初始条件（更保守的选择）
    x0 = [-1; -1; 10];  % 减小初始值
    
    % 时间序列
    nsteps = round(tmax/dt);
    t_span = (0:nsteps)*dt;
    
    % 初始化
    LE_history = zeros(n, nsteps+1);
    x = x0;
    Q = eye(n);
    
    % 缩短重正交化间隔
    orth_interval = 2;  % 更频繁的重正交化
    
    % 记录局部指数
    local_LEs = zeros(n, 1);
    steps_recorded = 0;
    
    try
        % 主循环
        for i = 1:nsteps
            % 计算系统演化
            [~, J] = NonlinearChaosSystem(t_span(i), x, params);
            
            % RK4积分（状态变量）
            k1 = NonlinearChaosSystem(t_span(i), x, params);
            k2 = NonlinearChaosSystem(t_span(i)+dt/2, x+dt*k1/2, params);
            k3 = NonlinearChaosSystem(t_span(i)+dt/2, x+dt*k2/2, params);
            k4 = NonlinearChaosSystem(t_span(i)+dt, x+dt*k3, params);
            x = x + dt*(k1 + 2*k2 + 2*k3 + k4)/6;
            
            % 切向空间演化
            Q_old = Q;
            for j = 1:n
                % 一阶方法用于切向量（简化计算以提高稳定性）
                Q(:,j) = Q(:,j) + dt * (J * Q(:,j));
            end
            
            % 定期重正交化
            if mod(i, orth_interval) == 0
                [Q, R] = qr(Q);  % QR分解
                
                % 计算局部指数
                for j = 1:n
                    local_LEs(j) = local_LEs(j) + log(abs(R(j,j)))/dt;
                end
                steps_recorded = steps_recorded + 1;
                
                % 计算当前的平均Lyapunov指数
                if steps_recorded > 0
                    LE_history(:,i+1) = local_LEs / steps_recorded;
                end
            else
                LE_history(:,i+1) = LE_history(:,i);
            end
            
            % 稳定性检查（使用更宽松的标准）
            if i > 100 && (any(abs(LE_history(:,i+1)) > 10) || ...
                         any(isnan(LE_history(:,i+1))))
                warning('检测到不稳定性，使用当前稳定结果');
                LE_history = LE_history(:,1:i);
                t_span = t_span(1:i);
                break;
            end
        end
        
        % 处理输出
        if size(LE_history,2) > 100
            % 使用后半段数据
            start_idx = floor(size(LE_history,2)/2);
            LE = LE_history(:,start_idx:end);
            t_span = t_span(start_idx:end);
            
            % 简单平滑处理
            window = min(20, floor(size(LE,2)/10));
            if window > 1
                for i = 1:n
                    LE(i,:) = movmean(LE(i,:), window);
                end
            end
        else
            LE = LE_history;
        end
        
    catch ME
        warning('计算过程出错: %s', ME.message);
        LE = zeros(n, 100);
        t_span = linspace(0, tmax, 100);
    end
    
    % 确保结果合理
    if size(LE,2) > 1
        % 对最终结果进行排序
        final_values = mean(LE(:,end-min(20,size(LE,2)-1):end), 2);
        [~, idx] = sort(final_values, 'descend');
        LE = LE(idx,:);
    end
end