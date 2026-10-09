function f = CECfn9(x)
    D = 30;
    split = floor(D/4);
    
    % Ackley
    part1 = -20*exp(-0.2*sqrt(sum(x(1:split).^2)/split))-exp(sum(cos(2*pi*x(1:split)))/split)+20+exp(1);
    % Rastrigin
    part2 = sum(x(split+1:2*split).^2 - 10*cos(2*pi*x(split+1:2*split))+10);
    % Griewank
    part3 = sum(x(2*split+1:3*split).^2)/4000 - prod(cos(x(2*split+1:3*split)./sqrt(1:split))) + 1;
    % Sphere
    part4 = sum(x(3*split+1:end).^2);
    
    f = part1 + part2 + part3 + part4 + 1100; % fmin = 1100
end