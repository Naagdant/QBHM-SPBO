function f = CECfn13(x)
    % CEC 2017 F15: Constrained Composition Function with 5 Components
    D = length(x);
    
    % Weighting factors
    w = [0.1, 0.25, 0.2, 0.15, 0.3];
    
    % Sphere
    part1 = sum(x.^2);
    
    % Rosenbrock
    part2 = sum(100 * (x(2:end) - x(1:end-1).^2).^2 + (x(1:end-1) - 1).^2);
    
    % Rastrigin
    part3 = sum(x.^2 - 10 * cos(2 * pi * x) + 10);
    
    % Griewank
    part4 = sum(x.^2) / 4000 - prod(cos(x ./ sqrt(1:D))) + 1;
    
    % Ackley
    part5 = -20 * exp(-0.2 * sqrt(sum(x.^2) / D)) - exp(sum(cos(2 * pi * x)) / D) + 20 + exp(1);
    
    % Weighted Sum of Components
    f = w(1) * part1 + w(2) * part2 + w(3) * part3 + w(4) * part4 + w(5) * part5 + 2000; % fmin = 2000
    
    % Constraint Handling: sum of squared values must not exceed a threshold
    constraint = sum(x.^2) - 200;
    if constraint > 0
        f = f + 1000 * constraint;
    end
end
