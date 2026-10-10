function [BestSol, BestFitness, Curve, PopHistRun] = bl_QPSO(PopSize, MaxIter, UB, LB, Dim, F_obj)
% QPSO (pure, no velocity) — compatible with test23.m
% Inputs/Outputs match: [BestSol, BestFitness, Curve, PopHistRun]

% ---- normalize bounds to 1xDim ----
if isscalar(UB), UB = repmat(UB,1,Dim); end
if isscalar(LB), LB = repmat(LB,1,Dim); end
range = UB - LB;

% ---- hyperparameters ----
c1 = 1.5; c2 = 1.5;            % cognitive/social weights (for local attractor)
alpha_max = 1.0;               % contraction factor start
alpha_min = 0.5;               % contraction factor end
alpha_nl  = 1.5;               % nonlinearity for alpha(t) decay

% ---- initialize population (same style as your bench) ----
X = zeros(PopSize, Dim);
for d = 1:Dim
    X(:,d) = LB(d) + (UB(d)-LB(d)) .* rand(PopSize,1);
end

% ---- evaluate initial & set archives ----
pBest = X;
pBest_fit = zeros(PopSize,1);
for i = 1:PopSize
    pBest_fit(i) = F_obj(X(i,:));
end
[BestFitness, ib] = min(pBest_fit);
BestSol = pBest(ib,:);

% ---- outputs ----
Curve = zeros(1, MaxIter);
PopHistRun = zeros(PopSize, Dim, MaxIter);

for iter = 1:MaxIter
    t = iter/MaxIter;
    % contraction factor alpha(t)
    alpha = alpha_min + (alpha_max - alpha_min) * (1 - t)^alpha_nl;

    % mBest: fitness-weighted mean of pBest (more stable than plain mean)
    w = 1 ./ (pBest_fit - min(pBest_fit) + 1e-12);
    w = w / sum(w);
    mBest = w' * pBest;              % 1xDim

    % local attractor p for each particle (canonical QPSO form)
    phi1 = rand(PopSize,Dim);
    phi2 = rand(PopSize,Dim);
    denom = (c1*phi1 + c2*phi2) + eps;
    P = (c1*phi1.*pBest + c2*phi2.*repmat(BestSol,PopSize,1)) ./ denom;

    % QPSO position update: X_new = P ± alpha*|mBest - X|*ln(1/U)
    U = rand(PopSize,Dim);
    sgn = sign(rand(PopSize,Dim) - 0.5); sgn(sgn==0) = 1;
    step = alpha .* abs(repmat(mBest,PopSize,1) - X) .* log(1./U);
    Xnew = P + sgn .* step;

    % bounds handling (reflection then clamp)
    Xnew = reflectBounds(Xnew, LB, UB);

    % evaluate & update archives
    for i = 1:PopSize
        f = F_obj(Xnew(i,:));
        if f < pBest_fit(i)
            pBest_fit(i) = f;
            pBest(i,:)   = Xnew(i,:);
        end
        if f < BestFitness
            BestFitness = f;
            BestSol     = Xnew(i,:);
        end
    end

    % advance swarm
    X = Xnew;

    % record
    Curve(iter) = BestFitness;
    PopHistRun(:,:,iter) = X;

    if mod(iter,50)==0 || iter==1
        div = mean(std(X,1,1) ./ (range + eps));
        fprintf('Iter %3d | Best: %.6e | alpha: %.3f | div: %.4f\n', ...
                iter, BestFitness, alpha, div);
    end
end
end

% ---------- helper ----------
function Xr = reflectBounds(X, LB, UB)
    LBm = repmat(LB,size(X,1),1);
    UBm = repmat(UB,size(X,1),1);
    Xr = X;
    below = Xr < LBm;  Xr(below) = LBm(below) + (LBm(below) - Xr(below));
    above = Xr > UBm;  Xr(above) = UBm(above) - (Xr(above) - UBm(above));
    Xr = min(max(Xr, LBm), UBm);
end

