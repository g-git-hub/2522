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
    
    % 系统方程 - 新的数学形式
    dx = -a*(x - y) + y*z;      % 基础耦合+非线性项
    dy = r*x - x*z - y;         % 能量输入-非线性抑制
    dz = x*y - b*z + c*x*y*z;   % 基础项+三维耦合
    
    % 状态限制
    x = min(max(x, -60), 60);
    y = min(max(y, -60), 60);
    z = min(max(z, 0), 80);    % z保持非负
    
    dX = [dx; dy; dz];
    
    if nargout > 1
        % Jacobian矩阵
        J = [
            -a,      a + z,     y;
            r - z,   -1,       -x;
            y + c*y*z, x + c*x*z, -b + c*x*y
        ];
    end
end