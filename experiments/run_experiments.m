% run_experiments  Reproduce all experiments of the paper (30 runs each).
% Every run is saved to results/<set>/<alg>_pXX_rYY.mat by cl_job (resumable).
% Add the repository folders to the MATLAB path first:
%   addpath(genpath(<repository root>))
% CEC 2022 needs the official code (cec22_test_func + input_data) from
%   https://github.com/P-N-Suganthan/2022-SO-BO  (not redistributed here).
root = pwd;
algs = {'QBHM','SPBO','PSO','QPSO','GWO','EO','TLBO','JAYA','GSA','FA'};
abl  = {'QBHM_orig','QBHM_noAbs','QBHM_noHawk'};   % ablation variants (Table 6b)
R = 30;

% ---- CECfn1-CECfn15, D = 30, 300 030 function evaluations ----
cfg = struct('root',root,'pop',30,'maxFE',300030,'D',30);
for a = [algs abl], for p = 1:15, for r = 1:R
    cl_job('CECfn', a{1}, p, r, cfg);
end, end, end

% ---- official CEC 2022, D = 10, 200 000 function evaluations ----
cfg = struct('root',root,'pop',30,'maxFE',200000,'D',10);
for a = [algs abl], for p = 1:12, for r = 1:R
    cl_job('CEC22', a{1}, p, r, cfg);
end, end, end

% ---- dispatch cases: {system, demand (MW), iterations of SPBO/QBHM-SPBO}, population 50 ----
cfg = struct('root',root,'pop',50,'maxFE',0,'D',0);
cfg.cases = {'ED3',2500,500; 'ED26',2000,500; 'ED26',2200,500; 'CEED10',2000,1000; ...
             'CEED6',150,500; 'CEED6',175,500; 'CEED6',200,500; 'CEED6',225,500; 'ED5',1800,500};
for a = [algs abl], for p = 1:9, for r = 1:R
    cl_job('DISP', a{1}, p, r, cfg);
end, end, end
