function chaos_system_analysis()
    % 系统参数
    params = [10.7644, 2.6261, 30.2043, 4.4504];
    
    % 参数设置
    r_range = linspace(25, 35, 500);
    dt = 0.01;
    tmax = 100;
    tspan = 0:dt:tmax;
    
    % 1. 分岔图
    fprintf('计算分岔图...\n');
    figure('Name', '分岔图', 'Position', [100, 100, 800, 400]);
    [r_points, x_points] = calculate_bifurcation(params, r_range, tspan);
    plot(r_points, x_points, '.b', 'MarkerSize', 0.3);
    xlabel('参数 r');
    ylabel('x 状态变量');
    title('分岔图');
    grid on;
    
    % 2. Lyapunov指数谱
    fprintf('计算Lyapunov指数谱...\n');
    figure('Name', 'Lyapunov指数谱', 'Position', [100, 550, 800, 400]);
    le_spectrum = calculate_le_spectrum(params, r_range);
    plot(r_range, le_spectrum, 'LineWidth', 1.5);
    xlabel('参数 r');
    ylabel('Lyapunov指数');
    title('Lyapunov指数谱');
    legend('\lambda_1', '\lambda_2', '\lambda_3', 'Location', 'best');
    grid on;
    
    % 3. Time Evolution
    fprintf('计算时间演化...\n');
    dt_final = 0.005;
    tspan_final = 0:dt_final:tmax;
    [t, x] = simulate_system(params, tspan_final);
    
    plot_time_evolution(t, x);
    
    fprintf('绘图完成!\n');
end

function [r_points, x_points] = calculate_bifurcation(params, r_range, tspan)
    x_points = [];
    r_points = [];
    
    parfor i = 1:length(r_range)
        temp_params = params;
        temp_params(3) = r_range(i);
        
        % 模拟系统
        [~, x] = simulate_system(temp_params, tspan);
        
        % 取后80%的数据点
        start_idx = floor(size(x,1)*0.2);
        x_data = x(start_idx:end,1);
        
        % 增加采样间隔以减少点数
        sample_interval = 8;
        sampled_points = x_data(1:sample_interval:end);
        
        % 获取极值点
        [max_peaks, ~] = findpeaks(x_data);
        [min_peaks, ~] = findpeaks(-x_data);
        min_peaks = -min_peaks;
        
        % 合并采样点和极值点
        x_all = [sampled_points; max_peaks; min_peaks];
        
        if ~isempty(x_all)
            temp_x = x_all;
            temp_r = ones(length(x_all),1)*r_range(i);
            
            x_points = [x_points; temp_x];
            r_points = [r_points; temp_r];
        end
    end
end

function le_spectrum = calculate_le_spectrum(params, r_range)
    le_spectrum = zeros(length(r_range), 3);
    
    parfor i = 1:length(r_range)
        temp_params = params;
        temp_params(3) = r_range(i);
        
        try
            % 计算李雅普诺夫指数
            tspan = 0:0.01:50;  % 用于计算指数的时间跨度
            [t, x] = simulate_system(temp_params, tspan);
            
            % 取后半段数据计算指数
            half_len = floor(length(t)/2);
            x = x(half_len:end,:);
            
            % 计算三个指数
            n = size(x,1);
            lambda = zeros(1,3);
            for j = 1:3
                dx = diff(x(:,j));
                dt = diff(t(half_len:end));
                lambda(j) = mean(log(abs(dx./dt)));
            end
            
            le_spectrum(i,:) = sort(lambda, 'descend');
        catch
            le_spectrum(i,:) = [0, 0, 0];
        end
    end
end

function [t, x] = simulate_system(params, tspan)
    x0 = [1; 1; 1];
    options = odeset('RelTol', 1e-4, 'AbsTol', 1e-7);
    [t, x] = ode45(@(t,x) NonlinearChaosSystem(t,x,params), tspan, x0, options);
end

function plot_time_evolution(t, x)
    % x分量
    figure('Name', 'Time Evolution - x', 'Position', [950, 100, 800, 250]);
    plot(t, x(:,1), 'LineWidth', 0.8);
    xlabel('时间'); ylabel('x');
    title('Time Evolution - x分量');
    grid on;
    
    % y分量
    figure('Name', 'Time Evolution - y', 'Position', [950, 400, 800, 250]);
    plot(t, x(:,2), 'LineWidth', 0.8);
    xlabel('时间'); ylabel('y');
    title('Time Evolution - y分量');
    grid on;
    
    % z分量
    figure('Name', 'Time Evolution - z', 'Position', [950, 700, 800, 250]);
    plot(t, x(:,3), 'LineWidth', 0.8);
    xlabel('时间'); ylabel('z');
    title('Time Evolution - z分量');
    grid on;
end

function [dX, J] = NonlinearChaosSystem(t, X, params)
    % 状态变量
    x = X(1);
    y = X(2);
    z = X(3);
    
    % 基本参数
    a = params(1);  % 耦合强度
    b = params(2);  % 阻尼系数
    r = params(3);  % 主控参数
    c = params(4);  % 非线性增益
    
    % 系统方程 - 创新修改
    dx = -a*(x - y) + y*z + 0.1*y^2;         % 增加y的二次项
    dy = r*x - x*z - y - 0.05*sin(z);        % 增加z的三角函数项
    dz = x*y - b*z + c*x*y*z + 0.1*x^2;      % 增加x的二次项
    
    % 状态限制 - 略微调整范围
    x = min(max(x, -65), 65);     % 扩大x范围
    y = min(max(y, -65), 65);     % 扩大y范围
    z = min(max(z, 0), 85);       % 扩大z范围，保持非负
    
    dX = [dx; dy; dz];
    
    if nargout > 1
        % 更新后的Jacobian矩阵
        J = [
            -a,      a + z + 0.2*y,     y;
            r - z,   -1,                -x - 0.05*cos(z);
            y + c*y*z + 0.2*x, x + c*x*z, -b + c*x*y
        ];
    end
end