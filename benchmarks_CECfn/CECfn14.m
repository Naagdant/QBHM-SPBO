function f = CECfn14(x)
    D = length(x);
    
    % Shifted transformation
    shift = 5;  
    x_shifted = x - shift;
    
    % Rastrigin function
    f_rastrigin = sum(x_shifted.^2 - 10 * cos(2 * pi * x_shifted) + 10);
    
    % Constraint: sum(x) ≤ 10 * D
    constraint = sum(x) - 10 * D;  
    penalty = max(0, constraint)^2 * 100; 
    
    % Final function value
    f = f_rastrigin + penalty + 1500;  % fmin = 1500
end
