function [BestSol, BestFitness, Curve, PopHistRun] = bl_GWO(PopSize, MaxIter, UB, LB, Dim, F_obj)
% GWO_wrapper  — standard GWO with strict FFE budget:
%   - PopSize points evaluated per iteration
%   - No extra offspring beyond the population
%   - Bounds enforced by projection
% Matches: [BestSol, BestFitness, Curve, PopHistRun]

    % ----- config -----
    doLocalSearch = false;   % set true if you want Nelder–Mead + box clamp
    % -------------------

    % vectorize bounds
    if isscalar(LB), LBv = LB*ones(1,Dim); else, LBv = LB(:).'; end
    if isscalar(UB), UBv = UB*ones(1,Dim); else, UBv = UB(:).'; end
    range = UBv - LBv;

    % init population
    X = LBv + rand(PopSize,Dim).*range;
    Fit = zeros(PopSize,1);
    for i=1:PopSize, Fit(i) = F_obj(X(i,:)); end

    % alpha, beta, delta wolves (best 3)
    [~,ord] = sort(Fit,'ascend');
    alphaX = X(ord(1),:); alphaF = Fit(ord(1));
    betaX  = X(ord(2),:); betaF  = Fit(ord(2));
    deltaX = X(ord(3),:); deltaF = Fit(ord(3));

    Curve = zeros(MaxIter,1);
    PopHistRun = cell(MaxIter,1);

    for it = 1:MaxIter
        a = 2 - 2*(it-1)/(MaxIter-1);  % linearly decreases from 2 to 0

        for i=1:PopSize
            r1 = rand(1,Dim); r2 = rand(1,Dim);
            A1 = 2*a*r1 - a; C1 = 2*r2;
            D_alpha = abs(C1.*alphaX - X(i,:));
            X1 = alphaX - A1.*D_alpha;

            r1 = rand(1,Dim); r2 = rand(1,Dim);
            A2 = 2*a*r1 - a; C2 = 2*r2;
            D_beta  = abs(C2.*betaX  - X(i,:));
            X2 = betaX  - A2.*D_beta;

            r1 = rand(1,Dim); r2 = rand(1,Dim);
            A3 = 2*a*r1 - a; C3 = 2*r2;
            D_delta = abs(C3.*deltaX - X(i,:));
            X3 = deltaX - A3.*D_delta;

            X(i,:) = (X1 + X2 + X3)/3;
        end

        % project to bounds & evaluate
        X = min(max(X,LBv),UBv);
        for i=1:PopSize, Fit(i) = F_obj(X(i,:)); end

        % update alpha/beta/delta
        [~,ord] = sort(Fit,'ascend');
        if Fit(ord(1)) < alphaF
            alphaX = X(ord(1),:); alphaF = Fit(ord(1));
        end
        if Fit(ord(2)) < betaF
            betaX  = X(ord(2),:); betaF  = Fit(ord(2));
        end
        if Fit(ord(3)) < deltaF
            deltaX = X(ord(3),:); deltaF = Fit(ord(3));
        end

        Curve(it) = alphaF;
        PopHistRun{it} = X;
    end

    % optional Nelder–Mead (keep box constraints inside objective)
    if doLocalSearch
        F_box = @(z) F_obj(min(max(z, LBv), UBv));
        options = optimset('Display','on','MaxIter',50,'TolX',1e-8,'TolFun',1e-8);
        [BestSol, BestFitness] = fminsearch(F_box, alphaX, options);
    else
        BestSol = alphaX; BestFitness = alphaF;
    end
end

