function [BestSol, BestFitness, Curve] = bl_EO(PopSize, MaxIter, LB, UB, Dim, F_obj)
% Wrapper for EO_noMemory (elitism OFF)

    % Ensure row vectors length Dim
    if isscalar(LB), lb = repmat(LB, 1, Dim); else, lb = LB(:)'; end
    if isscalar(UB), ub = repmat(UB, 1, Dim); else, ub = UB(:)'; end

    Run_no = 1;  % your main benchmark script loops over runs

    [Curve, ~, ~] = EO(PopSize, MaxIter, lb, ub, Dim, F_obj, Run_no);

    Curve = Curve(:)'; 
    BestFitness = Curve(end);

    % EO original code doesn't return best vector robustly; keep NaNs
    BestSol = nan(1, Dim);

    % Ensure curve length = MaxIter
    if isempty(Curve)
        Curve = nan(1, MaxIter);
    elseif numel(Curve) < MaxIter
        Curve = [Curve, repmat(Curve(end), 1, MaxIter-numel(Curve))];
    elseif numel(Curve) > MaxIter
        Curve = Curve(1:MaxIter);
    end
end
