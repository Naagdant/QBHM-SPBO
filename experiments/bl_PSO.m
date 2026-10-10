function [BestSol, BestFitness, Curve, PopHistRun] = bl_PSO(PopSize, MaxIter, UB, LB, Dim, F_obj)
% Plain PSO (gbest, Clerc constriction) with safe bound handling.

% --- bounds as row vectors ---
if isscalar(UB), UB = repmat(UB,1,Dim); end
if isscalar(LB), LB = repmat(LB,1,Dim); end
UB = UB(:)'; LB = LB(:)'; range = UB - LB;

% --- PSO params (Clerc constriction) ---
c1 = 2.05; c2 = 2.05;
phi = c1 + c2;
chi = 2 / (phi - 2 + sqrt(phi^2 - 4*phi));   % ~0.7298
Vmax = 0.20 * range;                          % per-dim velocity clamp

% --- init ---
rng('shuffle');
X = LB + rand(PopSize,Dim).*range;
V = zeros(PopSize,Dim);

fit = zeros(PopSize,1);
for i = 1:PopSize, fit(i) = F_obj(X(i,:)); end
pBest = X; pBest_fit = fit;
[BestFitness, gidx] = min(fit); BestSol = X(gidx,:);

% --- outputs ---
Curve = zeros(1,MaxIter);
PopHistRun = zeros(PopSize,Dim,MaxIter);

for iter = 1:MaxIter
    r1 = rand(PopSize,Dim);
    r2 = rand(PopSize,Dim);

    % velocity update (gbest)
    V = chi * ( V ...
        + c1.*r1 .* (pBest - X) ...
        + c2.*r2 .* (BestSol - X) );

    % clamp velocity per-dimension
    V = max(-Vmax, min(Vmax, V));

    % position update
    X = X + V;

    % ----- SAFE BOUND HANDLING (no shape-mismatch) -----
    % compute masks first (supports old MATLAB via bsxfun-like ops)
    hitLow  = bsxfun(@lt, X, LB);
    hitHigh = bsxfun(@gt, X, UB);

    % clamp positions with elementwise min/max using bsxfun (works everywhere)
    X = bsxfun(@min, bsxfun(@max, X, LB), UB);

    % zero velocity where we hit bounds (same-size logical mask)
    V(hitLow | hitHigh) = 0;

    % evaluate & update memories
    for i = 1:PopSize
        fi = F_obj(X(i,:));
        fit(i) = fi;
        if fi <= pBest_fit(i)
            pBest_fit(i) = fi;
            pBest(i,:)   = X(i,:);
            if fi < BestFitness
                BestFitness = fi;
                BestSol     = X(i,:);
            end
        end
    end

    Curve(iter) = BestFitness;
    PopHistRun(:,:,iter) = X;
end
end

