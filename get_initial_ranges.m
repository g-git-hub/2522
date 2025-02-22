function ranges = get_initial_ranges()
    ranges = struct();
    % 在基准值附近选择较小的扰动范围
    ranges.a = [-0.2, 0.2];      % 围绕 a=10
    ranges.b = [-0.1, 0.1];      % 围绕 b=8/3
    ranges.r = [-0.5, 0.5];      % 围绕 r=31.49
    ranges.c = [-0.2, 0.2];      % 围绕 c=4.99
end