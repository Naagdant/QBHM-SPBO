function f = CECfn15(x)
    D = length(x);
    
    % Shifted transformation
    shift = 3;  
    x_shifted = x - shift;
    
    % Griewank function
    f_griewank = sum(x_shifted.^2) / 4000 - prod(cos(x_shifted ./ sqrt(1:D))) + 1;
    
    % Constraint: Nonlinear equality constraint sum(x^2) = D
    constraint = abs(sum(x.^2) - D);
    penalty = constraint^2 * 1000; 
    
    % Final function value
    f = f_griewank + penalty + 2000;  % fmin = 2000
end
