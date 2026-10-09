function f = CECfn8(x)
    D = length(x);
    % Hybrid function combines Sphere, Schwefel, and Rosenbrock
    split = floor(D/3);
    part1 = sum(x(1:split).^2);                                   % Sphere
    part2 = sum(-x((1+floor(D/3)):2*floor(D/3)).*sin(sqrt(abs(x((1+floor(D/3)):2*floor(D/3)))))); % Schwefel
    part3 = sum(100*(x(2*floor(D/3)+1:D-1).^2 - x(2*floor(D/3)+2:D)).^2 + (x(2*floor(D/3)+1:D-1)-1).^2); % Rosenbrock
    f = part1 + part2 + part3 + 1000; % fmin = 1000
end