function f = CECfn3(x)
% CECfn3: Katsuura-type function (vectorised; identical to the loop version)
x = x(:)'; D = numel(x); j = (1:32)'; p = 2.^j;
t = sum(abs(p.*x - round(p.*x))./p, 1);
f = (10/D^2)*prod(1 + (1:D).*t) - 10/D^2 + 500;
end
