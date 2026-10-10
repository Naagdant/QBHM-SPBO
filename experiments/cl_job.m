function out = cl_job(setname, alg, prob, run, cfg)
% cl_job  One independent run, saved to disk (resumable).
%   setname : 'CECfn' | 'CEC22' | 'DISP'
%   alg     : 'SPBO','QBHM','QBHM_noHawk','QBHM_noAbs','QBHM_orig',
%             'PSO','GSA','JAYA','TLBO','FA','QPSO','GWO','EO'
%   prob    : function index (CECfn/CEC22) or case index (DISP)
%   cfg     : struct with pop, maxFE, D, cases (DISP), root
% Every algorithm receives exactly cfg.maxFE objective evaluations: the
% objective is wrapped by a counter that records the best value found within
% the budget and stops the run when the budget is exhausted.
root = cfg.root;
odir = fullfile(root,'results',setname); if ~exist(odir,'dir'), mkdir(odir); end
fn = fullfile(odir, sprintf('%s_p%02d_r%02d.mat', alg, prob, run));
if exist(fn,'file'), S = load(fn); out = S.out; return; end
cd(root);
seed = 2026 + run;  rng(seed,'twister');
pop = cfg.pop; maxFE = cfg.maxFE;

% ---------------- problem ----------------
det = [];
switch setname
  case 'CECfn'
    D = cfg.D; LB = -100*ones(1,D); UB = 100*ones(1,D);
    f0 = str2func(sprintf('CECfn%d',prob));
  case 'CEC22'
    D = cfg.D; LB = -100*ones(1,D); UB = 100*ones(1,D);
    f0 = @(x) cec22_test_func(x(:), prob);
  case 'DISP'
    c = cfg.cases(prob,:); pr = cl_dispatch(c{1}, c{2});
    D = pr.D; LB = pr.LB; UB = pr.UB; f0 = pr.fobj; det = pr;
    maxFE = pop*(c{3}+1);
end

% ---------------- budgeted objective ----------------
nFE = 0; bestF = inf; bestX = []; ck = round(linspace(1,maxFE,100)); conv = nan(1,100); ci = 1;
    function f = fb(x)
        if nFE >= maxFE, error('cl:budget','budget exhausted'); end
        x = min(max(x(:)',LB),UB);
        f = f0(x); nFE = nFE + 1;
        if f < bestF, bestF = f; bestX = x; end
        while ci <= 100 && nFE >= ck(ci), conv(ci) = bestF; ci = ci + 1; end
    end

% ---------------- algorithm ----------------
kEval = 1; if strcmp(alg,'TLBO'), kEval = 2; end
T = floor((maxFE - pop)/(pop*kEval));
div = []; h = @fb;
UBa = UB; LBa = LB; if all(UB==UB(1)) && all(LB==LB(1)), UBa = UB(1); LBa = LB(1); end
t0 = tic;
try
  switch alg
    case 'SPBO',        [~,~,~,info] = cl_SPBO(pop,T,UB,LB,D,h,maxFE); div = info.diversity;
    case 'QBHM',        [~,~,~,info] = QBHM_SPBO(pop,T,UB,LB,D,h,maxFE,struct()); div = info.diversity;
    case 'QBHM_noHawk', [~,~,~,info] = QBHM_SPBO(pop,T,UB,LB,D,h,maxFE,struct('hawking',false)); div = info.diversity;
    case 'QBHM_noAbs',  [~,~,~,info] = QBHM_SPBO(pop,T,UB,LB,D,h,maxFE,struct('absorb',false)); div = info.diversity;
    case 'QBHM_orig',   quiet(@cl_origQBHM, pop,T,UBa,LBa,D,h);
    case 'PSO',         quiet(@bl_PSO, pop,T,UBa,LBa,D,h);
    case 'GSA',         quiet(@bl_GSA, pop,T,UBa,LBa,D,h);
    case 'JAYA',        quiet(@bl_JAYA, pop,T,UBa,LBa,D,h);
    case 'TLBO',        quiet(@bl_TLBO, pop,T,UBa,LBa,D,h);
    case 'FA',          quiet(@bl_FA, pop,T,UBa,LBa,D,h);   % FA uses ~N^2/2 evals/iter: stopped by budget
    case 'QPSO',        quiet(@bl_QPSO, pop,T,UBa,LBa,D,h);
    case 'GWO',         quiet(@bl_GWO, pop,T,UBa,LBa,D,h);
    case 'EO',          quiet(@bl_EO, pop,T,LBa,UBa,D,h);
  end
catch ME
  if ~strcmp(ME.identifier,'cl:budget') && ~contains(ME.message,'budget exhausted'), rethrow(ME); end
end
out.time = toc(t0);
conv(ci:end) = bestF;
out.best = bestF; out.x = bestX; out.FE = nFE; out.curve = conv; out.alg = alg; out.prob = prob; out.run = run; out.seed = seed;
if ~isempty(div), out.div = div(round(linspace(1,find(~isnan(div),1,'last'),100))); end
if ~isempty(det), out.det = det.detail(bestX); out.best = out.det.total; end
save(fn,'out');
end

function quiet(fun, varargin)
evalc('fun(varargin{:});');
end
