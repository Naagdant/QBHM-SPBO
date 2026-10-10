function [Best_sol, Best_fit, curve, pop_hist_QBHM_SPBO] = ...
         cl_origQBHM(pop, max_iter, UB, LB, D, fobj)
% QBHM-SPBO: Quantum Black Hole Model–based SPBO
%
% USAGE:
%   [Best_sol, Best_fit, curve, pop_hist] = ...
%       cl_origQBHM(pop, max_iter, UB, LB, D, fobj);
%
% Inputs:
%   pop       – population size
%   max_iter  – max number of iterations
%   UB, LB    – upper/lower bounds (scalar or 1×D)
%   D         – problem dimension
%   fobj      – handle to objective function (fobj(x))
%
% Outputs:
%   Best_sol    – best solution found
%   Best_fit    – its fitness
%   curve       – 1×max_iter convergence history
%   pop_hist    – pop×D×max_iter array of populations

% ensure LB and UB are 1×D
if isscalar(LB), LB = repmat(LB,1,D); end
if isscalar(UB), UB = repmat(UB,1,D); end

% 1) Initialization
X   = LB + rand(pop, D) .* (UB - LB);
fit = zeros(pop,1);
for i = 1:pop
    fit(i) = fobj(X(i,:));
end

[Best_fit, idx] = min(fit);
Best_sol        = X(idx,:);

% preallocate
curve              = zeros(1, max_iter);
pop_hist_QBHM_SPBO = [];

% 2) Main loop
for iter = 1:max_iter
    mean_X = mean(X,1);
    radius = Best_fit / (sum(fit) + eps);

    X_new = zeros(pop, D);  % new population buffer

    for i = 1:pop
        rand_student = randi(pop);
        eta          = rand();  % quantum threshold

        for d = 1:D
            if fit(i) == Best_fit
                % --- Update best student ---
                if eta <= radius
                    % Quantum Black Hole absorption
                    alpha = adaptive_alpha(iter, max_iter);
                    X_new(i,d) = Best_sol(d) + ...
                        alpha * abs(Best_sol(d) - X(rand_student,d)) * log(1./rand());
                else
                    % Regular SPBO best‐student update
                    m = randi([1,2]);
                    X_new(i,d) = Best_sol(d) + ...
                        (-1)^m * rand() * (Best_sol(d) - X(rand_student,d));
                end
            else
                % --- SPBO for non‐best students ---
                if rand() < rand()
                    if rand()>0.5
                        % Good student (Eq.2a)
                        X_new(i,d) = Best_sol(d) + ...
                            rand() * (Best_sol(d) - X(i,d));
                    else
                        % Good student (Eq.2b)
                        X_new(i,d) = X(i,d) + ...
                            rand() * (Best_sol(d) - X(i,d)) + ...
                            rand() * (X(i,d) - mean_X(d));
                    end
                else
                    if rand()>0.5
                        % Average student (Eq.3)
                        X_new(i,d) = X(i,d) + ...
                            rand() * (mean_X(d) - X(i,d));
                    else
                        % Random improvement (Eq.4)
                        X_new(i,d) = LB(d) + rand() * (UB(d) - LB(d));
                    end
                end
            end

            % --- Boundary clamp (per‐dimension) ---
            X_new(i,d) = min( max( X_new(i,d), LB(d) ), UB(d) );
        end

        % evaluate and greedy select
        new_fit = fobj(X_new(i,:));
        if new_fit < fit(i)
            X(i,:) = X_new(i,:);
            fit(i) = new_fit;
        end
    end

    % update global best
    [current_best, idx] = min(fit);
    if current_best < Best_fit
        Best_fit = current_best;
        Best_sol = X(idx,:);
    end

    % record
    curve(iter)                  = Best_fit;
    

    % display progress

end

end  % function QBHM_SPBO

%% Sub‐function: adaptive expansion–contraction factor
function alpha = adaptive_alpha(t, Tmax)
    alpha_max = 1.5;
    alpha_min = 0.1;
    alpha     = alpha_max - ((alpha_max - alpha_min) * (t / Tmax));
end
