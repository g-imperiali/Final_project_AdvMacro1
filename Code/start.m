clear;
clc;
close all
addpath("../CompEcon2017 silicon/CEtools")
intmeth      = 'spline'; % interpolation: 'linear' or 'spline'
optset('bisect', 'tol', 1e-12);

% Assigned Parameters
p.theta   = 2;
p.eta   = 10;
p.theta_h   = 2.02;   % from reduced form estimates
p.eta_h   = 10;

p.xi      = 5; 

p.S         = 50000;                            % number of sectors
p.N         = 5;                                 % number of firms in sector

rng(0); 
z         = (1-rand(p.S,p.N)).^(-1/p.xi);      % generate S by N matrix of draws from Pareto(xi)

% define helpful exponents
exp1 = p.eta_h*(1-p.eta)/(p.eta+p.eta_h);
exp2 = -p.eta*(1+p.eta_h)/(p.eta+p.eta_h);

% find mu(x) function. Easier to find inverse
mumin     = p.eta/(p.eta - 1);
mumax     = p.theta/(p.theta - 1);
dmin     = p.eta_h/(p.eta_h - 1); % d for delta
dmax     = p.theta_h/(p.theta_h - 1);

mux       = (nodeunif(100000, 0, 1).^3)*(mumax - mumin) + mumin;    % more nodes in low markup region

x         = ((1/p.theta - 1/p.eta).*mux.^(exp1)).^(-1).*(1 - 1./mux - 1/p.eta);
p.xmin = min(x);
p.xmax = max(x);
p.Mu      = griddedInterpolant({x}, mux, intmeth);  % interpolate mu(x)

dg       = (nodeunif(100000, 0, 1).^3)*(dmax - dmin) + dmin;    % more nodes in low markup region
g        = ((1/p.theta_h - 1/p.eta_h).*dg.^(exp2)).^(-1).*(dg -1 - 1/p.eta_h);
p.gmin = min(g);
p.gmax = max(g);
p.D      = griddedInterpolant({g}, dg, intmeth);  % interpolate delta(x)

% Iteration
d_old = ones(p.S,p.N);

for itt = 1 : 100

Dsold   = mean(z, 2); % guess for DS1 and DS2
Ds1 = solve_broyden('solve_DS1', Dsold, 1e-12*ones(p.S,1), 10*ones(p.S,1), p,z, d_old, exp1);

x = (z.^(-1/p.eta_h-1).*d_old./Ds1).^(exp1);
x = min(max(x, p.xmin), p.xmax);
mu_old = p.Mu(x);
%errmu = solve_mu(mu, p,x);

Ds2 = solve_broyden('solve_DS2', Dsold, 1e-6*ones(p.S,1), 10*ones(p.S,1), p,z, mu_old,exp2);
g = (z.^(1/p.eta-1).*mu_old./Ds2).^(exp2);
g = min(max(g, p.gmin), p.gmax);
d_new = p.D(g);

% recalculate mu
Ds1 = solve_broyden('solve_DS1', Ds1, 1e-12*ones(p.S,1), 10*ones(p.S,1), p,z, d_new,exp1);
err   = solve_DS1(Ds1, p, z,d_new,exp1) ;% check if the equation that determines D is satisfied (S by 1)
x = (z.^(-1/p.eta_h-1).*d_new./Ds1).^(exp1);
x = min(max(x, p.xmin), p.xmax);
mu_new = p.Mu(x);

parold    = [d_old; mu_old];
parnew    = [d_new; mu_new];

fprintf('\n')
fprintf('iteration: %d\n', itt)
fprintf('error in industry equilibrium   = %9.3e \n',  norm(err));
fprintf('error in firm markups           = %9.3e \n', norm(parnew - parold));

fprintf('\n')

step     = 1/2;

d_old     = d_new + step*(d_new - d_old);
mu_old    = mu_new + step*(mu_new - mu_old);

if norm(parnew - parold) < 1e-7

    break
end


end

% Aggregation

