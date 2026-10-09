function f = CECfn7(x)
    D = 30;
    f = 0;
    for i = 1:D-1
        temp = x(i)^2 + x(i+1)^2;
        f = f + 0.5 + (sin(sqrt(temp))^2 - 0.5)/(1 + 0.001*temp)^2;
    end
    f = f + 900; % fmin = 900
end