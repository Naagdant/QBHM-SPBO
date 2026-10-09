function f = CECfn5(x)
    D = 30;
    sum_x = sum(x); sum_x2 = sum(x.^2);
    term1 = sqrt(abs(sum_x2^2 - sum_x^2));
    term2 = (0.5*sum_x2 + sum_x)/D;
    f = term1 + term2 + 700; % fmin = 700
end