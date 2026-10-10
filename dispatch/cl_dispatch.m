function prob = cl_dispatch(name, PD)
% cl_dispatch  Clean dispatch test systems for the QBHM-SPBO revision.
%   prob = cl_dispatch('ED3',2500) | ('ED26',2000|2200) | ('CEED10',2000) | ('CEED6',150..225)
% Returns:
%   prob.fobj(x)  objective for the optimiser (x = outputs of units 2..n)
%   prob.LB, prob.UB, prob.D
%   prob.detail(x) struct: P (full vector), fuel, emission, loss, total, violation
% Power balance: unit 1 is the slack unit, P1 = PD + PL(P) - sum(P2..Pn),
% solved by fixed-point iteration when losses are present. If P1 leaves its
% limits it is clamped, the balance error |sum(P) - PD - PL| is reported as
% the violation and penalised with rho*violation^2 (rho = 1e6), identical for
% every algorithm.

rho = 1e6;
switch upper(name)
  case 'ED3'   % 3-unit cubic-cost ED, no losses
    a = [749.55 1285 1531]'; b = [6.95 7.05 6.531]';
    c = [9.68e-4 7.375e-4 1.04e-3]'; d = [1.27e-7 6.453e-8 9.98e-8]';
    Pmin = [320 300 275]'; Pmax = [800 1200 1100]';
    fuelF = @(P) a + b.*P + c.*P.^2 + d.*P.^3;
    emisF = @(P) zeros(size(P)); hpen = 0; B = [];
  case 'ED5'   % 5-unit cubic-cost ED (units 1,2,3 of ED3 plus copies of units 1 and 2), no losses
    a = [749.55 1285 1531 749.55 1285]'; b = [6.95 7.05 6.531 6.95 7.05]';
    c = [9.68e-4 7.375e-4 1.04e-3 9.68e-4 7.375e-4]'; d = [1.27e-7 6.453e-8 9.98e-8 1.27e-7 6.453e-8]';
    Pmin = [320 300 275 320 300]'; Pmax = [800 1200 1100 800 1200]';
    fuelF = @(P) a + b.*P + c.*P.^2 + d.*P.^3;
    emisF = @(P) zeros(size(P)); hpen = 0; B = [];
  case 'ED26'  % 26-unit cubic-cost ED, no losses
    a = [24.38 24.41 24.63 24.76 24.88 117.75 118.10 118.45 118.82 81.13 81.29 81.46 81.62 ...
         217.89 218.33 218.77 142.73 143.02 143.31 143.59 259.13 259.64 260.17 177.05 310.00 311.91]';
    b = [25.54 25.67 25.80 25.93 26.06 37.55 37.66 37.77 37.88 13.32 13.35 13.38 13.40 ...
         18.00 18.09 18.20 10.69 10.71 10.73 10.75 23.00 23.10 23.20 10.86 7.49 7.50]';
    c = [0.025 0.026 0.028 0.028 0.028 0.011 0.012 0.013 0.014 0.008 0.008 0.009 0.009 ...
         0.006 0.006 0.005 0.004 0.004 0.004 0.004 0.002 0.002 0.002 0.001 0.001 0.001]';
    d = [5.08e-9 -1.01e-8 1.01e-8 -5.08e-9 -5.72e-16 8.31e-8 8.56e-8 8.15e-8 8.29e-8 ...
         -5.80e-10 -5.47e-10 -5.49e-10 -5.50e-10 1.25e-18 -1.19e-18 2.44e-18 1.11e-10 ...
         1.03e-10 1.03e-10 1.03e-10 1.07e-10 1.04e-10 1.00e-10 -4.42e-19 -1.10e-19 -3.55e-20]';
    Pmin = [2.4*ones(5,1); 4*ones(4,1); 15.2*ones(4,1); 25*ones(3,1); 54.25*ones(4,1); 68.95*ones(3,1); 140; 100; 100];
    Pmax = [12*ones(5,1); 20*ones(4,1); 76*ones(4,1); 100*ones(3,1); 155*ones(4,1); 197*ones(3,1); 350; 400; 400];
    % slack must be a large unit: reorder so unit 24 is the slack
    fuelF = @(P) a + b.*P + c.*P.^2 + d.*P.^3;
    emisF = @(P) zeros(size(P)); hpen = 0; B = [];
  case 'CEED10' % 10-unit CEED, valve-point fuel, exponential emission, B-loss
    G = [0.12951 40.5407 1000.403 33 0.0174 10 55; 0.10908 39.5804 950.606 25 0.0178 20 80;
         0.12511 36.5104 900.705 32 0.0162 47 120; 0.12111 39.5104 800.705 30 0.0168 20 130;
         0.15247 38.5390 756.799 30 0.0148 50 160; 0.10587 46.1592 451.325 20 0.0163 70 240;
         0.03546 38.3055 1243.531 20 0.0152 60 300; 0.02803 40.3965 1049.998 30 0.0128 70 340;
         0.02111 36.3278 1658.569 60 0.0136 135 470; 0.01799 38.2704 1356.659 40 0.0141 150 470];
    E = [0.04702 -3.9864 360.0012 0.25475 0.01234; 0.04652 -3.9524 350.0056 0.25475 0.01234;
         0.04652 -3.9023 330.0056 0.25163 0.01215; 0.04652 -3.9023 330.0056 0.25163 0.01215;
         0.00420 0.3277 13.8593 0.24970 0.01200;  0.00420 0.3277 13.8593 0.24970 0.01200;
         0.00680 -0.5455 40.2669 0.24800 0.01290; 0.00680 -0.5455 40.2669 0.24990 0.01203;
         0.00460 -0.5112 42.8955 0.25470 0.01234; 0.00460 -0.5112 42.8955 0.25470 0.01234];
    B = 1e-4*[0.49 0.14 0.15 0.15 0.16 0.17 0.17 0.18 0.19 0.20; 0.14 0.45 0.16 0.16 0.17 0.15 0.15 0.16 0.18 0.18;
              0.15 0.16 0.39 0.10 0.12 0.12 0.14 0.14 0.16 0.16; 0.15 0.16 0.10 0.40 0.14 0.10 0.11 0.12 0.14 0.15;
              0.16 0.17 0.12 0.14 0.35 0.11 0.13 0.13 0.15 0.16; 0.17 0.15 0.12 0.10 0.11 0.36 0.12 0.12 0.14 0.15;
              0.17 0.15 0.14 0.11 0.12 0.12 0.38 0.16 0.16 0.18; 0.18 0.16 0.14 0.12 0.13 0.12 0.16 0.40 0.15 0.16;
              0.19 0.18 0.16 0.14 0.15 0.14 0.16 0.15 0.42 0.19; 0.20 0.18 0.16 0.15 0.16 0.15 0.18 0.16 0.19 0.44];
    Pmin = G(:,6); Pmax = G(:,7);
    fuelF = @(P) G(:,1).*P.^2 + G(:,2).*P + G(:,3) + abs(G(:,4).*sin(G(:,5).*(Pmin - P)));
    emisF = @(P) E(:,1).*P.^2 + E(:,2).*P + E(:,3) + E(:,4).*exp(E(:,5).*P);
    hpen  = fuelF(Pmax)./emisF(Pmax);          % max-max price penalty factor per unit ($/lb)
  case 'CEED6'  % IEEE 30-bus 6-unit CEED, cubic fuel and CO2/SO2/NOx, no losses
    A  = [-136.0 14.50 0.0920 0.0010; -3.50 22.00 0.0250 0.0004; -81.00 23.00 0.0750 0.0006;
          -14.50 13.50 0.1000 0.0002; -9.75 11.50 0.1200 0.0013; 75.60 12.50 0.0840 0.0004];
    SO = [-90.00 17.00 0.150 0.0005; -30.50 12.00 0.055 0.0014; -80.00 10.00 0.035 0.0010;
          -34.50 23.50 0.070 0.0020; -19.75 21.50 0.120 0.0013; -25.60 22.50 0.080 0.0021];
    NO = [-26.00 18.50 0.0520 0.0012; -35.00 12.00 0.0450 0.0004; -15.00 13.00 0.0500 0.0016;
          -74.50 17.50 0.0700 0.0012; -89.05 8.50 0.0400 0.0003; -75.00 15.50 0.0240 0.0014];
    CO = [-16.00 14.00 0.0920 0.0015; -93.50 12.50 0.0250 0.0014; -85.00 13.50 0.0550 0.0016;
          -24.50 13.50 0.0100 0.0012; -59.00 21.00 0.0400 0.0023; -70.00 22.00 0.0800 0.0014];
    hSO = [1.0851 1.061 2.1051 0.5975 0.6771 0.6470]';
    hNO = [0.9407 1.4961 1.3869 0.8308 2.1705 1.0930]';
    hCO = [0.7822 1.1895 1.4355 1.1332 0.7455 0.7158]';
    Pmin = [50 20 15 10 10 12]'; Pmax = [200 80 50 50 50 40]';
    cub = @(K,P) K(:,1) + K(:,2).*P + K(:,3).*P.^2 + K(:,4).*P.^3;
    fuelF = @(P) cub(A,P);
    emisF = @(P) [cub(CO,P) cub(SO,P) cub(NO,P)];   % kg/h, columns CO2 SO2 NOx
    hpen  = [hCO hSO hNO]; B = [];
  otherwise
    error('unknown system');
end

n = numel(Pmin);
if strcmpi(name,'ED26'), sl = 24; elseif strcmpi(name,'ED5'), sl = 3; else, sl = 1; end   % slack unit index
oth = setdiff(1:n, sl);

prob.name = name; prob.PD = PD; prob.n = n; prob.slack = sl;
prob.LB = Pmin(oth)'; prob.UB = Pmax(oth)'; prob.D = n-1;
prob.Pmin = Pmin; prob.Pmax = Pmax; prob.hpen = hpen;
prob.fobj   = @(x) objective(x);
prob.detail = @(x) detail(x);
prob.evalP  = @(P) evalFull(P(:));

    function P = buildP(x)
        P = zeros(n,1); P(oth) = x(:);
        P(sl) = PD - sum(P(oth));
        if ~isempty(B)
            for it = 1:20
                P(sl) = PD + P'*B*P - sum(P(oth));
            end
        end
        P(sl) = min(max(P(sl), Pmin(sl)), Pmax(sl));
    end
    function s = evalFull(P)
        s.P = P;
        s.loss = 0; if ~isempty(B), s.loss = P'*B*P; end
        s.fuel = sum(fuelF(P));
        Em = emisF(P);
        s.emission = sum(Em,1);
        s.emisCost = sum(hpen.*Em,1);
        s.total = s.fuel + sum(s.emisCost);
        s.violation = abs(sum(P) - PD - s.loss);
        s.limitViol = sum(max(0,Pmin-P) + max(0,P-Pmax));
    end
    function f = objective(x)
        s = evalFull(buildP(x));
        f = s.total + rho*(s.violation^2 + s.limitViol^2);
    end
    function s = detail(x)
        s = evalFull(buildP(x));
    end
end
