function f = CECfn4(x)
    D = 30;
    sum_x = sum(x); sum_x2 = sum(x.^2);
    term1 = abs(sum_x2 - D)^(1/4);
    term2 = (0.5*sum_x2 + sum_x)/D;
    f = term1 + term2 + 600; % fmin = 600
end