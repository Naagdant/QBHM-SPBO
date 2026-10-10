function [Best_sol, Best_fit, curve, info] = cl_SPBO(pop, max_iter, UB, LB, D, fobj, maxFE)
% cl_SPBO  Student Psychology-Based Optimization (Das et al., 2020).
% Population-synchronous implementation of Eqs. (1)-(4); identical update
% rules to the original SPBO.m. Same interface/budget as cl_QBHMSPBO for a
% fair comparison (one function evaluation per student per iteration).
%
% maxFE (optional) : stop after this many function evaluations (default pop*max_iter + pop)

if nargin < 7 || isempty(maxFE), maxFE = pop*max_iter + pop; end
if isscalar(LB), LB = repmat(LB,1,D); end
if isscalar(UB), UB = repmat(UB,1,D); end
LBm = repmat(LB,pop,1);  UBm = repmat(UB,pop,1);

X   = LB + rand(pop,D).*(UB-LB);
fit = evalPop(fobj, X);  FE = pop;
[Best_fit, ib] = min(fit);  Best_sol = X(ib,:);
curve = nan(1,max_iter);  divc = nan(1,max_iter);
tStart = tic;

for t = 1:max_iter
    if FE + pop > maxFE, break; end
    meanX = mean(X,1);
    isBest = (fit == Best_fit);  isBest(find(isBest,1)+1:end) = false; % one best student
    Xn = X;

    % ---- best student (Eq. 1) ----
    ib = find(isBest);
    if ~isempty(ib)
        k  = randi(pop,1,D);
        sg = (-1).^randi([1 2],1,D);
        Xn(ib,:) = Best_sol + sg.*rand(1,D).*(Best_sol - X(sub2ind([pop D],k,1:D)));
    end

    % ---- other students (Eqs. 2a,2b,3,4), decided per dimension ----
    nb  = ~isBest;  Xo = X(nb,:);  m = size(Xo,1);
    good = rand(m,D) < rand(m,D);
    sub  = rand(m,D) > 0.5;
    B  = repmat(Best_sol,m,1);  M = repmat(meanX,m,1);
    e2a = B + rand(m,D).*(B - Xo);
    e2b = Xo + rand(m,D).*(B - Xo) + rand(m,D).*(Xo - M);
    e3  = Xo + rand(m,D).*(M - Xo);
    e4  = LBm(nb,:) + rand(m,D).*(UBm(nb,:) - LBm(nb,:));
    Y = e4;
    Y( good &  sub) = e2a( good &  sub);
    Y( good & ~sub) = e2b( good & ~sub);
    Y(~good &  sub) = e3(~good &  sub);
    Xn(nb,:) = Y;

    Xn = min(max(Xn,LBm),UBm);
    fn = evalPop(fobj, Xn);  FE = FE + pop;
    imp = fn < fit;  X(imp,:) = Xn(imp,:);  fit(imp) = fn(imp);

    [cb, ib] = min(fit);
    if cb < Best_fit, Best_fit = cb; Best_sol = X(ib,:); end
    curve(t) = Best_fit;
    divc(t)  = popDiversity(X, LB, UB);
end
info.FE = FE;  info.time = toc(tStart);  info.diversity = divc;
end

function f = evalPop(fobj, X)
f = zeros(size(X,1),1);
for i = 1:size(X,1), f(i) = fobj(X(i,:)); end
end

function dv = popDiversity(X, LB, UB)
% mean Euclidean distance to centroid, normalised by search-space diagonal
dv = mean(sqrt(sum((X - mean(X,1)).^2,2))) / norm(UB - LB);
end
