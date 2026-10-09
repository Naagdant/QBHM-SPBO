function f = CECfn11(x)
    D = 30;
    split = floor(D/3);

    % Ackley
    part1 = -20*exp(-0.2*sqrt(mean(x(1:split).^2))) - exp(mean(cos(2*pi*x(1:split)))) + 20 + exp(1);
    
    % Griewank
    part2 = sum(x(split+1:2*split).^2)/4000 - prod(cos(x(split+1:2*split)./sqrt(1:split))) + 1;
    
    % Schwefel
    part3 = sum(-x(2*split+1:D).*sin(sqrt(abs(x(2*split+1:D)))));
    
    f = part1 + part2 + part3 + 1400; % fmin = 1400
end
