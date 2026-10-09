function f = CECfn6(x)
    D = 30;
    z = x + 1 + 420.9687462275036;
    f = 0;
    for i = 1:D-1
        f = f + 100*(z(i)^2 - z(i+1))^2 + (z(i) - 1)^2;
    end
    f = f - 800; % Shift to fmin = 800
end