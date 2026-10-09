% demo.m  Minimal example of QBHM-SPBO on two standard test functions.
% Usage: run this script in MATLAB (R2018b or later); no toolbox is required.
clear; clc; rng(1,'twister');

D     = 30;            % problem dimension
pop   = 30;            % number of students
maxFE = 10000*D;       % function-evaluation budget
T     = floor((maxFE - pop)/pop);   % iterations (one evaluation per student per iteration)
LB = -100; UB = 100;

problems = { 'Sphere',    @(x) sum(x.^2);
             'Rastrigin', @(x) sum(x.^2 - 10*cos(2*pi*x) + 10) };

figure('Color','w');
for k = 1:size(problems,1)
    [x, fbest, curve, info] = QBHM_SPBO(pop, T, UB, LB, D, problems{k,2}, maxFE);
    fprintf('%-10s best f = %.4e   (%d evaluations, %.1f s)\n', problems{k,1}, fbest, info.FE, info.time);
    subplot(1,2,k); semilogy(curve,'LineWidth',1.2); grid on;
    xlabel('Iteration'); ylabel('Best objective value'); title(problems{k,1});
end

% Options (all optional), e.g. ablation variants used in the paper:
%   QBHM_SPBO(pop,T,UB,LB,D,f,maxFE,struct('hawking',false))  % without Hawking re-emission
%   QBHM_SPBO(pop,T,UB,LB,D,f,maxFE,struct('absorb',false))   % without quantum absorption
