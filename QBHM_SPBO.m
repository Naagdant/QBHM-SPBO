function [Best_sol, Best_fit, curve, info] = QBHM_SPBO(pop, max_iter, UB, LB, D, fobj, maxFE, opts)
% QBHM_SPBO  Quantum-inspired black hole model based Student Psychology-Based Optimization (QBHM-SPBO).
%
% The best student acts as a quantum-inspired black hole. A student whose
% normalised fitness lies inside the event horizon (Schwarzschild radius
% R_t) is "absorbed": it is resampled from a quantum delta-potential well
% centred between itself and the black hole (exploitation). Students
% outside the horizon keep the original SPBO Eqs. (2)-(4) (exploration).
% Students that stagnate for L iterations "evaporate" (Hawking radiation):
% they are re-emitted by an opposition jump about the black hole.
%
% Normalised fitness  s_i = (f_i - f_best)/(f_worst - f_best + eps)  in [0,1]
% is invariant to shifting/scaling of f, so the absorption condition is
% valid for zero, negative and zero-centred objectives.
%
% Absorption condition (single, used everywhere):   s_i <= R_t
%   R_t = R_min + (R_max - R_min) * t/T        (horizon grows -> exploitation late)
% Quantum move (per dimension):
%   phi ~ U(0,1),  p_d = phi*Xbest_d + (1-phi)*X_i,d
%   X_i,d <- p_d +/- beta_t * |Xmean_d - X_i,d| * ln(1/u),  u ~ U(0,1)
%   beta_t = beta_max - (beta_max - beta_min) * t/T   (contraction-expansion)
%
% opts (all optional, for ablation):
%   .absorb    (true)  quantum absorption operator on/off
%   .hawking   (true)  Hawking re-emission on/off
%   .Rmin .Rmax        horizon range          (0.10, 0.50)
%   .bmax .bmin        beta range             (1.00, 0.50)
%   .L                 stagnation limit       (max(20, round(T/50)))

if nargin < 7 || isempty(maxFE), maxFE = pop*max_iter + pop; end
if nargin < 8, opts = struct(); end
o = defaults(opts, max_iter);
if isscalar(LB), LB = repmat(LB,1,D); end
if isscalar(UB), UB = repmat(UB,1,D); end
LBm = repmat(LB,pop,1);  UBm = repmat(UB,pop,1);

X   = LB + rand(pop,D).*(UB-LB);
fit = evalPop(fobj, X);  FE = pop;
[Best_fit, ib] = min(fit);  Best_sol = X(ib,:);
stag  = zeros(pop,1);
curve = nan(1,max_iter);  divc = nan(1,max_iter);
nAbs = zeros(1,max_iter); nHaw = zeros(1,max_iter);
tStart = tic;

for t = 1:max_iter
    if FE + pop > maxFE, break; end
    tau  = t/max_iter;
    Rt   = o.Rmin + (o.Rmax - o.Rmin)*tau;
    beta = o.bmax - (o.bmax - o.bmin)*tau;
    meanX = mean(X,1);
    fw = max(fit);
    s  = (fit - Best_fit) ./ (fw - Best_fit + eps);        % normalised fitness

    isBest = (fit == Best_fit);  isBest(find(isBest,1)+1:end) = false;
    inH = (s <= Rt) & ~isBest;                             % absorption condition
    if ~o.absorb, inH(:) = false; end
    Xn = X;

    % ---- best student: SPBO Eq. (1) ----
    ib = find(isBest);
    if ~isempty(ib)
        k  = randi(pop,1,D);
        sg = (-1).^randi([1 2],1,D);
        Xn(ib,:) = Best_sol + sg.*rand(1,D).*(Best_sol - X(sub2ind([pop D],k,1:D)));
    end

    % ---- absorbed students: quantum delta-well collapse ----
    ia = find(inH);  na = numel(ia);
    if na > 0
        Xa  = X(ia,:);
        phi = rand(na,D);
        P   = phi.*repmat(Best_sol,na,1) + (1-phi).*Xa;
        sg  = sign(rand(na,D) - 0.5);
        Xn(ia,:) = P + sg.*beta.*abs(repmat(meanX,na,1) - Xa).*log(1./rand(na,D));
    end

    % ---- remaining students: SPBO Eqs. (2a)-(4) ----
    nb = ~isBest & ~inH;  Xo = X(nb,:);  m = size(Xo,1);
    if m > 0
        good = rand(m,D) < rand(m,D);  sub = rand(m,D) > 0.5;
        B = repmat(Best_sol,m,1);  M = repmat(meanX,m,1);
        e2a = B + rand(m,D).*(B - Xo);
        e2b = Xo + rand(m,D).*(B - Xo) + rand(m,D).*(Xo - M);
        e3  = Xo + rand(m,D).*(M - Xo);
        e4  = LBm(nb,:) + rand(m,D).*(UBm(nb,:) - LBm(nb,:));
        Y = e4;
        Y( good &  sub) = e2a( good &  sub);
        Y( good & ~sub) = e2b( good & ~sub);
        Y(~good &  sub) = e3(~good &  sub);
        Xn(nb,:) = Y;
    end

    % ---- Hawking re-emission of stagnant students (not the best) ----
    ih = [];
    if o.hawking
        ih = find(stag >= o.L & ~isBest);
        if ~isempty(ih)
            nh = numel(ih);
            % opposition point about the black hole, randomly scaled
            Xn(ih,:) = repmat(Best_sol,nh,1) + rand(nh,D).*(repmat(Best_sol,nh,1) - X(ih,:));
            stag(ih) = 0;
        end
    end

    Xn = min(max(Xn,LBm),UBm);
    fn = evalPop(fobj, Xn);  FE = FE + pop;
    imp = fn < fit;
    X(imp,:) = Xn(imp,:);  fit(imp) = fn(imp);
    stag(imp) = 0;  stag(~imp) = stag(~imp) + 1;

    [cb, ib] = min(fit);
    if cb < Best_fit, Best_fit = cb; Best_sol = X(ib,:); end
    curve(t) = Best_fit;  divc(t) = popDiversity(X, LB, UB);
    nAbs(t) = na;  nHaw(t) = numel(ih);
end
info.FE = FE;  info.time = toc(tStart);  info.diversity = divc;
info.nAbsorbed = nAbs;  info.nHawking = nHaw;  info.opts = o;
end

function o = defaults(o, T)
d = struct('absorb',true,'hawking',true,'Rmin',0.10,'Rmax',0.50, ...
           'bmax',1.00,'bmin',0.50,'L',max(20,round(T/50)));
fn = fieldnames(d);
for k = 1:numel(fn), if ~isfield(o,fn{k}), o.(fn{k}) = d.(fn{k}); end, end
end

function f = evalPop(fobj, X)
f = zeros(size(X,1),1);
for i = 1:size(X,1), f(i) = fobj(X(i,:)); end
end

function dv = popDiversity(X, LB, UB)
dv = mean(sqrt(sum((X - mean(X,1)).^2,2))) / norm(UB - LB);
end
