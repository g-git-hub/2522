function [LE, T] = compute_lyapunov_wolf(params)
    % 优化的初始化参数
    dt = 0.005;          % 增大时间步长
    tspan = 0:dt:30;     % 缩短总计算时间
    n = 3;               % 系统维度
    GS_step = 10;        % 减少正交化频率
    
    % 初始条件
    x0 = [1; 1; 1];
    d0 = 1e-7;          % 略微增大初始扰动
    v = eye(n) * d0;
    
    % 预分配内存
    num_steps = length(tspan);
    LE = zeros(num_steps, n);
    T = tspan;
    
    % 状态初始化
    x = x0;
    sum_log = zeros(n, 1);
    steps = 0;
    
    % 优化的ODE选项
    options = odeset('RelTol', 1e-6, 'AbsTol', 1e-6);
    
    try
        for i = 1:num_steps
            % 单步演化
            [~, Y] = ode45(@(t,y) fast_system_evolution(t, y, params), ...
                [0, dt], [x; reshape(v, n*n, 1)], options);
            
            % 提取新状态
            x_new = Y(end, 1:n)';
            v_new = reshape(Y(end, n+1:end), n, n);
            
            % 定期正交化
            if mod(i, GS_step) == 0
                [Q, R] = qr(v_new);
                sum_log = sum_log + log(abs(diag(R))/d0);
                v_new = Q * diag(ones(n,1) * d0);
                steps = steps + 1;
                
                % 计算当前LE
                if steps > 0
                    LE(i,:) = sum_log / (steps * GS_step * dt);
                end
            else
                LE(i,:) = LE(max(1,i-1),:);  % 复制上一步的值
            end
            
            % 更新状态
            x = x_new;
            v = v_new;
        end
        
        % 结果验证和排序
        if any(isnan(LE(end,:))) || any(isinf(LE(end,:)))
            error('Invalid LE');
        end
        LE = sort(LE, 2, 'descend');
        
    catch
        LE = zeros(num_steps, n) - 1;
    end
end

function dy = fast_system_evolution(t, y, params)
    n = 3;
    x = y(1:n);
    v = reshape(y(n+1:end), n, n);
    
    % 快速计算系统动力学和Jacobian
    [dx, J] = quick_system_calc(x, params);
    
    % 线性化方程
    dv = J * v;
    
    % 组合结果
    dy = [dx; reshape(dv, n*n, 1)];
end

function [dx, J] = quick_system_calc(x, params)
    % 直接计算，避免函数调用开销
    a = params(1);
    b = params(2);
    r = params(3);
    c = params(4);
    
    % 系统方程
    dx = [-a*(x(1) - x(2)) + x(2)*x(3);
          r*x(1) - x(1)*x(3) - x(2);
          x(1)*x(2) - b*x(3) + c*x(1)*x(2)*x(3)];
    
    % Jacobian矩阵
    J = [-a,        a + x(3),    x(2);
         r - x(3),  -1,         -x(1);
         x(2) + c*x(2)*x(3), x(1) + c*x(1)*x(3), -b + c*x(1)*x(2)];
end