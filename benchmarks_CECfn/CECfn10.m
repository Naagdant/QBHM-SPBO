function f = CECfn10(x)
    D = 30;
    split = floor(D/5);

    % Sphere
    part1 = sum(x(1:split).^2);
    
    % Rosenbrock
    part2 = sum(100*(x(split+1:2*split-1).^2 - x(split+2:2*split)).^2 + (x(split+1:2*split-1)-1).^2);

    % Ackley
    part3 = -20*exp(-0.2*sqrt(mean(x(2*split+1:3*split).^2))) - exp(mean(cos(2*pi*x(2*split+1:3*split)))) + 20 + exp(1);
    
    % Rastrigin
    part4 = sum(x(3*split+1:4*split).^2 - 10*cos(2*pi*x(3*split+1:4*split)) + 10);
    
    % Schwefel
    part5 = sum(-x(4*split+1:end).*sin(sqrt(abs(x(4*split+1:end)))));
    
    f = part1 + part2 + part3 + part4 + part5 + 1300; % fmin = 1300
end
