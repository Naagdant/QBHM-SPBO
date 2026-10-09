function f = CECfn12(x)
    % CEC 2017 F11: Constrained Hybrid Composition Function
    D = length(x);
    w = [0.3, 0.2, 0.5];  % Weighting factors for sub-functions
    
    % Rastrigin
    part1 = sum(x.^2 - 10 * cos(2 * pi * x) + 10);
    
    % Griewank
    part2 = sum(x.^2) / 4000 - prod(cos(x ./ sqrt(1:D))) + 1;
    
    % Ackley
    part3 = -20 * exp(-0.2 * sqrt(sum(x.^2) / D)) - exp(sum(cos(2 * pi * x)) / D) + 20 + exp(1);
    
    % Weighted Sum of Sub-Functions
    f = w(1) * part1 + w(2) * part2 + w(3) * part3 + 1500; % fmin = 1500
    
    % Constrained Handling (adds penalty if sum of x is greater than a threshold)
    constraint = sum(x) - 5;
    if constraint > 0
        f = f + 1000 * constraint;  % Large penalty
    end
end
