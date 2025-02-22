function [r_values, x_peaks] = bifurcation_diagram(param_range, param_index, num_points)
    % 输入参数检查
    if nargin < 3
        num_points = 1000;
    end
    
    % 初始化
    r_values = linspace(param_range(1), param_range(2), num_points);
    x_peaks = cell(num_points, 1);
    
    % 基准参数
    base_params = [0, 0, 0, 0];
    
    % 计算每个参数值的峰值
    parfor i = 1:num_points
        % 设置当前参数
        params = base_params;
        params(param_index) = r_values(i);
        
        % 模拟系统
        tspan = [0 500];  % 足够长的时间
        x0 = [-5 + 10*rand(1);
              -5 + 10*rand(1);
              15 + 10*rand(1)];
        
        options = odeset('RelTol', 1e-8, 'AbsTol', 1e-8);
        [t, X] = ode45(@(t,x) NonlinearChaosSystem(t,x,params), tspan, x0, options);
        
        % 去除暂态
        start_idx = round(length(t)/2);
        X = X(start_idx:end, 1);  % 只取x分量
        
        % 找出局部极大值
        peaks = findpeaks(X);
        x_peaks{i} = peaks;
    end
end