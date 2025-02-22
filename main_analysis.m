% 主程序：改进的混沌系统分析
clear all;
close all;

try
    % 1. 参数优化
    [best_params, best_le] = optimize_parameters(30);
    
    if ~all(isinf(best_le))
        fprintf('\n最终优化结果:\n');
        fprintf('最佳参数:\n');
        fprintf('a = %.4f\n', best_params(1));
        fprintf('b = %.4f\n', best_params(2));
        fprintf('r = %.4f\n', best_params(3));
        fprintf('c = %.4f\n', best_params(4));
        
        % 2. 计算完整的Lyapunov谱
        [LE, t] = compute_lyapunov_wolf(best_params);
        
        % 绘制结果
        figure('Name', 'Lyapunov Spectrum');
        plot(t, LE, 'LineWidth', 1.5);
        xlabel('Time');
        ylabel('Lyapunov Exponents');
        grid on;
        title('Lyapunov Spectrum Evolution');
        legend('\lambda_1', '\lambda_2', '\lambda_3');
        
        % 3. 相空间分析
        tspan = [0 30];
        x0 = [-1; -1; 10];
        options = odeset('RelTol', 1e-8, 'AbsTol', 1e-8);
        [t, X] = ode45(@(t,x) NonlinearChaosSystem(t,x,best_params), ...
                       tspan, x0, options);
        
        figure('Name', 'Phase Space');
        plot3(X(:,1), X(:,2), X(:,3));
        grid on;
        xlabel('x');
        ylabel('y');
        zlabel('z');
        title('Optimized Chaotic Attractor');
        
        % 4. 时间序列分析
        figure('Name', 'Time Series');
        subplot(3,1,1);
        plot(t, X(:,1));
        ylabel('x(t)');
        title('State Variables Evolution');
        grid on;
        
        subplot(3,1,2);
        plot(t, X(:,2));
        ylabel('y(t)');
        grid on;
        
        subplot(3,1,3);
        plot(t, X(:,3));
        xlabel('Time');
        ylabel('z(t)');
        grid on;
    else
        error('参数优化失败，未找到有效解');
    end
    
catch ME
    fprintf('错误: %s\n', ME.message);
end