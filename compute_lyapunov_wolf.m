function [LE, t_span] = compute_lyapunov_wolf(params)
    % 系统参数
    n = 3;          % 系统维数
    tstart = 0;     % 起始时间
    tend = 30;      % 减少总时间以提高稳定性
    stept = 0.05;   % 减小步长提高精度
    
    % 初始条件
    x0 = [-1; -1; 10];
    
    % 扩展系统初始化
    y = zeros(n*(n+1), 1);
    y(1:n) = x0;
    for i = 1:n
        y(n*i + i) = 1.0;  % 初始化单位矩阵
    end
    
    % 初始化存储
    nsteps = round((tend-tstart)/stept);
    t_span = zeros(nsteps, 1);
    LE = zeros(n, nsteps);
    cum = zeros(n, 1);
    
    try
        % 主循环
        t = tstart;
        for iter = 1:nsteps
            % 积分一步
            [~, Y] = ode45(@(t,y) extended_system(t,y,params), [t t+stept], y);
            
            if isempty(Y) || any(isnan(Y(end,:)))
                error('数值积分失败');
            end
            
            y = Y(end,:)';
            t = t + stept;
            
            % 提取切向向量
            Q = reshape(y(n+1:end), n, n);
            
            % QR分解
            [Q, R] = qr(Q);
            
            % 更新累积和
            for i = 1:n
                if abs(R(i,i)) > 1e-10  % 避免取log(0)
                    cum(i) = cum(i) + log(abs(R(i,i)));
                end
            end
            
            % 计算当前指数
            LE(:,iter) = cum / (t - tstart);
            t_span(iter) = t;
            
            % 更新变分方程初值
            y(n+1:end) = Q(:);
        end
        
    catch ME
        error('计算过程出错: %s', ME.message);
    end
    
    % 确保输出维度正确
    LE = LE';
end

function dy = extended_system(t, y, params)
    n = 3;
    x = y(1:n);
    Y = reshape(y(n+1:end), n, n);
    
    % 计算原系统和Jacobian
    [dx, J] = NonlinearChaosSystem(t, x, params);
    
    % 变分方程
    dY = J * Y;
    
    % 返回扩展系统
    dy = [dx; dY(:)];
end