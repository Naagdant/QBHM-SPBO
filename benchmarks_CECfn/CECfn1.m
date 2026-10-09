function f = CECfn1(x)
% CECfn1: Weierstrass-type function (vectorised; identical to the loop version)
D = numel(x); a = 0.5; b = 3; k = (0:20)'; ak = a.^k; bk = b.^k;
sum1 = sum(sum(ak .* cos(2*pi*bk .* (x(:)' + 0.5))));
sum2 = sum(ak .* cos(pi*bk));
f = sum1 - D*sum2 + 300;
end
