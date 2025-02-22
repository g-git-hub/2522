function plot_system_analysis()
    % 系统参数
    params = [10.7644, 2.6261, 30.2043, 4.4504, 0.1];  % 增加了新参数d=0.1
    
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
            [LE, ~] = compute_lyapunov_wolf(temp_params);
            le_spectrum(i,:) = mean(LE(end-5:end,:), 1);
        catch
            le_spectrum(i,:) = [0, 0, 0];
        end
    end
end

function [t, x] = simulate_system(params, tspan)
    x0 = [1; 1; 1];
    options = odeset('RelTol', 1e-4, 'AbsTol', 1e-7);
    [t, x] = ode45(@(t,x) system_equations(t,x,params), tspan, x0, options);
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

function dx = system_equations(t, x, params)
    a = params(1);
    b = params(2);
    r = params(3);
    c = params(4);
    d = params(5);  % 新增参数
    
    dx = zeros(3,1);
    % 修改后的方程组，增加了新的非线性耦合项
    dx(1) = -a*(x(1) - x(2)) + x(2)*x(3) + d*x(1)*x(2);  % 增加x1x2项
    dx(2) = r*x(1) - x(1)*x(3) - x(2) - d*sin(x(3));     % 增加sin(x3)项
    dx(3) = x(1)*x(2) - b*x(3) + c*x(1)*x(2)*x(3) + d*x(2)^2;  % 增加x2^2项
end