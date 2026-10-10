function [Best_sol, Best_fit, curve] = bl_JAYA(pop, max_iter, UB, LB, D, fobj)
% JAYA Algorithm
% Inputs:
% pop      - Population size
% max_iter - Maximum number of iterations
% UB, LB   - Upper and lower bounds
% D        - Dimensionality of the problem
% fobj     - Objective function handle
%
% Outputs:
% Best_sol - Best solution found
% Best_fit - Best fitness value
% curve    - Convergence curve

% Initialization
X = LB + rand(pop, D).*(UB - LB);
fit = zeros(pop,1);

for i=1:pop
    fit(i) = fobj(X(i,:));
end

[Best_fit, idx] = min(fit);
Best_sol = X(idx,:);
curve = zeros(1, max_iter);

for iter = 1:max_iter

    [best_fit_iter, best_idx] = min(fit);
    [worst_fit_iter, worst_idx] = max(fit);

    Best_X = X(best_idx,:);
    Worst_X = X(worst_idx,:);

    for i = 1:pop
        for d = 1:D
            r1 = rand();
            r2 = rand();

            % Jaya Algorithm update
            X_new(i,d) = X(i,d) + r1*(Best_X(d) - abs(X(i,d))) - r2*(Worst_X(d) - abs(X(i,d)));

            % Boundary check
            X_new(i,d) = min(max(X_new(i,d), LB(min(d,numel(LB)))), UB(min(d,numel(UB))));  % per-dimension bounds (identical for scalar bounds)
        end

        % Fitness evaluation
        new_fit = fobj(X_new(i,:));

        % Greedy selection
        if new_fit < fit(i)
            X(i,:) = X_new(i,:);
            fit(i) = new_fit;
        end
    end

    % Update best solution
    [current_best, idx] = min(fit);
    if current_best < Best_fit
        Best_fit = current_best;
        Best_sol = X(idx,:);
    end

    curve(iter) = Best_fit;

    % Optional progress display
    if mod(iter,100)==0 || iter==1
        fprintf('Iter %d: Best Fitness = %.6f\n', iter, Best_fit);
    end
end

end
