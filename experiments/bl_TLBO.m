function [Best_sol, Best_fit, curve] = bl_TLBO(pop, max_iter, UB, LB, D, fobj)
% TLBO - Teaching-Learning-Based Optimization
% Input parameters:
% pop       - Population size
% max_iter  - Maximum number of iterations
% UB, LB    - Upper and lower bounds (scalar)
% D         - Number of dimensions
% fobj      - Handle to the objective function
%
% Returns:
% Best_sol  - Best solution found
% Best_fit  - Best fitness value
% curve     - Convergence curve (best fitness per iteration)

% Initialize population
X = LB + rand(pop, D).*(UB - LB);
fit = zeros(pop,1);
for i=1:pop
    fit(i) = fobj(X(i,:));
end

[Best_fit, idx] = min(fit);
Best_sol = X(idx,:);
curve = zeros(1,max_iter);

for iter = 1:max_iter
    mean_X = mean(X,1);
    [~, Teacher_idx] = min(fit);
    Teacher = X(Teacher_idx,:);

    % Teaching Phase
    for i = 1:pop
        TF = randi([1,2]); % Teaching factor randomly selected as 1 or 2
        X_new = X(i,:) + rand(1,D).*(Teacher - TF*mean_X);

        % Boundary constraints
        X_new = min(max(X_new, LB), UB);

        new_fit = fobj(X_new);
        if new_fit < fit(i)
            X(i,:) = X_new;
            fit(i) = new_fit;
        end
    end

    % Learning Phase
    for i = 1:pop
        partner_idx = randi(pop);
        while partner_idx == i
            partner_idx = randi(pop);
        end

        X_new = X(i,:);
        if fit(i) < fit(partner_idx)
            X_new = X(i,:) + rand(1,D).*(X(i,:) - X(partner_idx,:));
        else
            X_new = X(i,:) + rand(1,D).*(X(partner_idx,:) - X(i,:));
        end

        % Boundary constraints
        X_new = min(max(X_new, LB), UB);

        new_fit = fobj(X_new);
        if new_fit < fit(i)
            X(i,:) = X_new;
            fit(i) = new_fit;
        end
    end
end
    % Update global best
