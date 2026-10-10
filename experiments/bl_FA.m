function [Best_sol, Best_fit, curve] = bl_FA(pop, max_iter, UB, LB, D, fobj)
% FA - Firefly Algorithm
% Inputs:
% pop       - Population size
% max_iter  - Maximum number of iterations
% UB, LB    - Upper and lower bounds (scalar)
% D         - Number of dimensions (variables)
% fobj      - Objective function handle
%
% Outputs:
% Best_sol  - Best solution found
% Best_fit  - Best fitness value
% curve     - Convergence curve

% FA parameters
alpha = 0.2;       % Randomness parameter
beta0 = 1;         % Attraction coefficient base value
gamma = 1;         % Light absorption coefficient

% Initialize fireflies randomly
X = LB + rand(pop, D).*(UB - LB);
fit = zeros(pop,1);
for i = 1:pop
    fit(i) = fobj(X(i,:));
end

[Best_fit, idx] = min(fit);
Best_sol = X(idx,:);
curve = zeros(1,max_iter);

for iter = 1:max_iter
    for i = 1:pop
        for j = 1:pop
            if fit(j) < fit(i)
                r = norm(X(i,:) - X(j,:)); % Euclidean distance
                beta = beta0 * exp(-gamma * r.^2); % Attraction factor
                % Update firefly position
                X(i,:) = X(i,:) + beta*(X(j,:) - X(i,:)) + alpha*(rand(1,D)-0.5).*(UB - LB);
                % Boundary check
                X(i,:) = min(max(X(i,:), LB), UB);
                % Evaluate fitness
                fit(i) = fobj(X(i,:));
            end
        end
    end

    % Update best solution
    [current_best, idx] = min(fit);
    if current_best < Best_fit
        Best_fit = current_best;
        Best_sol = X(idx,:);
    end

    curve(iter) = Best_fit;

    % Display progress
    if mod(iter,100)==0 || iter==1
        fprintf('Iteration %d: Best Fitness = %.6f\n', iter, Best_fit);
    end
end

end
