function [dX, J] = NonlinearChaosSystem(t, X, params)
    % 状态变量
    x = X(1);
    y = X(2);
    z = X(3);
    
    % 基准参数值
    a = 10;        % 耦合强度
    b = 8/3;       % 衰减速度
    r = 31.49;     % 对流强度参数
    c = 4.99;      % 非线性强度
    
    % 添加小扰动
    a = a + params(1);
    b = b + params(2);
    r = r + params(3);
    c = c + params(4);
    
    % 系统方程
    dx = -a*(x - y) + y*z;
    dy = -x*z + r*x - y;
    dz = x*y - b*z + c*x*y*z;
    
    dX = [dx; dy; dz];
    
    % Jacobian矩阵
    if nargout > 1
        J = [
            -a,      a+z,    y;
            r-z,    -1,     -x;
            y+c*y*z, x+c*x*z, -b+c*x*y
        ];
    end
end