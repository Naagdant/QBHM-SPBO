function f = CECfn2(x)
% CECfn2: Adjusted Schwefel’s Function
D = numel(x);
sum2 = 0;
for i = 1:D
    z = x(i) + 4.209687462275036e+002;
    if abs(z) <= 500
        sum_part = z*sin(sqrt(abs(z)));
    elseif z > 500
        sum_part = (500 - mod(z,500))*sin(sqrt(abs(500 - mod(z,500)))) ...
                  - ((z - 500)^2)/(10000*D);
    else
        sum_part = (mod(abs(z),500) - 500)*sin(sqrt(abs(mod(abs(z),500) - 500)));
    end
    sum2 = sum2 + sum_part; % Accumulate sum_part across all dimensions
end
f = 418.9829*D - sum2 + 400; % Use the accumulated sum
end