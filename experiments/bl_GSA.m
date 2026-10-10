function [Best_sol, Best_fit, curve] = bl_GSA(pop, max_iter, UB, LB, D, fobj)
% GSA - Gravitational Search Algorithm
%
% This algorithm is inspired by Newtonian gravity and mass interactions.
% Agents (solutions) attract each other proportional to their masses,
% which are derived from their fitness values.
%
% Inputs:
%   pop      - Population size (number of agents)
%   max_iter - Maximum number of iterations
%   UB, LB   - Upper and lower bounds (assumed to be scalars)
%   D        - Problem dimensionality (number of variables)
%   fobj     - Objective function handle (minimization)
%
% Outputs:
%   Best_sol - Best solution found
%   Best_fit - Best fitness value
%   curve    - Convergence curve (best fitness at each iteration)
%
% Initialization: random positions and zero velocities
X = LB + rand(pop, D) .* (UB - LB);
v = zeros(pop, D);  % Initial velocities are set to zero

% Evaluate initial fitness of all agents
fit = zeros(pop, 1);
for i = 1:pop
    fit(i) = fobj(X(i, :));
end
[Best_fit, bestIdx] = min(fit);
Best_sol = X(bestIdx, :);
curve = zeros(1, max_iter);

% GSA specific parameters
G0 = 100;         % Initial gravitational constant
alpha = 20;       % Gravitational constant decay coefficient
epsilon = 1e-6;   % Small constant to avoid division by zero

for iter = 1:max_iter
    % Update gravitational constant with exponential decay
    G = G0 * exp(-alpha * iter / max_iter);
    
    % Determine best and worst fitness for mass calculation
    best = min(fit);
    worst = max(fit);
    
    % Calculate mass for each agent based on its fitness (minimization)
    m = zeros(pop, 1);
    for i = 1:pop
        m(i) = (worst - fit(i)) / (worst - best + epsilon);
    end
    % Normalize masses so that they sum to one
    M = m / (sum(m) + epsilon);
    
    % Compute gravitational forces for each agent
    force = zeros(pop, D);
    for i = 1:pop
        for j = 1:pop
            if i ~= j
                % Compute Euclidean distance between agents i and j
                dist = norm(X(i,:) - X(j,:)) + epsilon;
                % Randomized factor to introduce stochasticity
                rand_factor = rand();
                % Accumulate force component from agent j on agent i
                force(i,:) = force(i,:) + rand_factor * G * (M(i) * M(j) / dist) * (X(j,:) - X(i,:));
            end
        end
    end
    
    % Update acceleration, velocity, and positions of agents
    a = zeros(pop, D);
    for i = 1:pop
        % Acceleration is the force divided by the agent's mass
        a(i,:) = force(i,:) / (M(i) + epsilon);
        % Update velocity: introduce inertia with a random factor
        v(i,:) = rand() .* v(i,:) + a(i,:);
        % Update position based on new velocity
        X(i,:) = X(i,:) + v(i,:);
        % Enforce boundary constraints
        X(i,:) = min(max(X(i,:), LB), UB);
    end
    
    % Evaluate new fitness values for all agents
    for i = 1:pop
        fit(i) = fobj(X(i, :));
    end
    
    % Update global best if a better solution is found
    [current_best, idx] = min(fit);
    if current_best < Best_fit
        Best_fit = current_best;
        Best_sol = X(idx, :);
    end
    
    % Record the best fitness in the convergence curve
    curve(iter) = Best_fit;
    
    % Optional: display progress every 100 iterations or at the first iteration
    if mod(iter,100)==0 || iter==1
        fprintf('Iteration %d: Best Fitness = %.6f\n', iter, Best_fit);
    end
end
end

