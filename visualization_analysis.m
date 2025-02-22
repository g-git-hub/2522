function visualization_analysis(params)
    % 生成时间序列数据
    tspan = [0 100];
    x0 = [-1; -1; 10];
    options = odeset('RelTol', 1e-8, 'AbsTol', 1e-8);
    [t, X] = ode45(@(t,x) NonlinearChaosSystem(t,x,params), tspan, x0, options);
    
    % 1. 三维混沌吸引子
    figure('Name', '3D Chaotic Attractor');
    plot3(X(:,1), X(:,2), X(:,3), 'b-', 'LineWidth', 0.5);
    grid on;
    xlabel('x');
    ylabel('y');
    zlabel('z');
    title('3D Chaotic Attractor');
    view(45, 30);
    
    % 2. x-y平面投影
    figure('Name', 'X-Y Plane');
    plot(X(:,1), X(:,2), 'b-', 'LineWidth', 0.5);
    grid on;
    xlabel('x');
    ylabel('y');
    title('X-Y Plane Projection');
    
    % 3. x-z平面投影
    figure('Name', 'X-Z Plane');
    plot(X(:,1), X(:,3), 'b-', 'LineWidth', 0.5);
    grid on;
    xlabel('x');
    ylabel('z');
    title('X-Z Plane Projection');
    
    % 4. y-z平面投影
    figure('Name', 'Y-Z Plane');
    plot(X(:,2), X(:,3), 'b-', 'LineWidth', 0.5);
    grid on;
    xlabel('y');
    ylabel('z');
    title('Y-Z Plane Projection');
    
    % 5. 时间序列
    figure('Name', 'Time Series');
    subplot(3,1,1);
    plot(t, X(:,1), 'b-');
    grid on;
    ylabel('x(t)');
    title('Time Evolution');
    
    subplot(3,1,2);
    plot(t, X(:,2), 'b-');
    grid on;
    ylabel('y(t)');
    
    subplot(3,1,3);
    plot(t, X(:,3), 'b-');
    grid on;
    xlabel('Time');
    ylabel('z(t)');
end