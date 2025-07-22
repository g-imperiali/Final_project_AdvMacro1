% =========================================================================
% Final project: "Two-sided market power and misallocation" 
% Authors      : Grace Chuan, Guglielmo Imperiali 
% Date         : July 18, 2025
% =========================================================================

clear; clc; close all

% Add path to CompEcon library
addpath(genpath('CompEcon'));

% % % Export command window results
% diary('ps11_results.txt');
% diary on;

% Create output folder for figures
figFolder = fullfile('../', 'figures/');
if ~exist(figFolder, 'dir')
    mkdir(figFolder);
end

% Set plot graphic parameters
set(groot, 'DefaultAxesLineWidth', 1.5);
set(groot, 'DefaultLineLineWidth', 3);
set(groot, 'DefaultAxesTickLabelInterpreter','latex'); 
set(groot, 'DefaultLegendInterpreter','latex');
set(groot, 'DefaultAxesFontSize',22);

printr       = 0;                  % if 1, report iteration errors

% Assigned Parameters

p.sigma   = 2;                     % CRRA
p.phi     = 1;                     % Frisch

p.theta   = 2;                     % CES on goods across sectors
p.eta     = 10;                    % CES on goods within sectors

p.thetaL  = 3;                     % CES on labor across sectors
p.etaL    = 5;                     % CES on labor within sectors

p.xi      = 5;                     % Pareto shape param

S         = 500000;                % number of sectors
NN        = [1, 2, 5, 10, 100];    % numbers of firms in sector

rng(0); 

N         = NN(3);                 % 5 firms per sector

fprintf('\n')
fprintf('Number of firms in each sector  = %d \n',  N);
fprintf('\n')

u         = rand(S, N);          % S-by-N uniformly distributed draws in [0,1]
z         = (1 - u).^(-1/p.xi);  % generate S-by-N matrix of draws from Pareto(xi)

%% Equilibrium within each sector

% Bounds on mu and grid
mumin     = p.eta/(p.eta - 1);     % when market share is 0
mumax     = p.theta/(p.theta - 1); % when market share is 1

mux       = (nodeunif(100000, 0, 1).^3)*(mumax - mumin) + mumin;            % more nodes in low markup region

% Bounds on delta and grid
deltamin  = (p.etaL + 1)/p.etaL;     % when wage bill share is 0
deltamax  = (p.thetaL + 1)/p.thetaL; % when wage bill share is 1

deltax    = (nodeunif(100000, 0, 1).^3)*(deltamax - deltamin) + deltamin;

% Guess f and g
exp1      = p.etaL*(1-p.eta)/(p.etaL + p.eta);  % exponent on mu in markup equation
exp2      = -p.eta*(1+p.etaL)/(p.etaL + p.eta); % exponent on delta in markdown equation

f         = ((1/p.theta - 1/p.eta) .* mux.^exp1).^(-1) .* ...
            (1 - 1./mux - 1/p.eta);

g         = abs( ((1/p.thetaL - 1/p.etaL) .* deltax.^exp2).^(-1) .* ...
            (deltax - 1 - 1/p.etaL) );          % Take absolute value to ensure all values are positive

p.Mu      = griddedInterpolant(log(f), mux, 'spline', 'nearest');
p.Delta   = griddedInterpolant(log(g), deltax, 'spline', 'nearest');  

% find D1(s) and D2(s) 
A1        = z.^(-(p.etaL + 1)/p.etaL);
A2        = z.^(-(p.eta  - 1)/p.eta);

optset('bisect', 'tol', 1e-18);

tic;

for i = 1:1000

    if ~exist('delta', 'var')
        delta = ones(S,N)/N; % starting guess for delta is uniform 1/N within each sector
    end
    
    deltaold = delta;
    
    % Given deltaold, compute mu 
    D1       = bisect('findD1', 1e-3, 1000, p, A1, exp1, deltaold); 
    [~, mu, ~, ~] = findD1(D1, p, A1, exp1, deltaold);

    % Update delta based on computed mu
    D2       = bisect('findD2', 1e-3, 1000, p, A2, exp2, mu); 
    [~, delta, ~, ~] = findD2(D2, p, A2, exp2, mu);
    
    if norm(delta - deltaold)/norm(deltaold) < 1e-15; break; end
    
    if i == 1000; disp('Warning: Algorithm did not converge'); end

end

time = toc;
fprintf('Computation time                = %.4f \n',  time);

% Compute final quantities with solutions from iterative algorithm
D1     = sum((mu .* delta .* A1).^exp1, 2) .^(1/exp1);
f      = (delta .* A1 ./ D1).^exp1; 

mu     = p.Mu(log(f));
omega  = f .* mu.^exp1;
err_mu = 1 - 1./mu - 1/p.eta - (1/p.theta - 1/p.eta) .* omega;

D2     = sum((mu .* delta .* A2).^exp2, 2) .^(1/exp2);
g      = (mu .* A2 ./ D2).^exp2;

delta  = p.Delta(log(g));
alpha  = g .* delta.^exp2;
err_delta = delta - 1 - 1/p.etaL - (1/p.thetaL - 1/p.etaL) .* alpha;

fprintf('\n')
fprintf('CONVERGENCE CHECK:                      \n');
fprintf('error in industry equilibrium   = %9.3e \n',  norm(sum(omega, 2) - 1));
fprintf('error in firm markups           = %9.3e \n',  norm(err_mu(:)));
fprintf('\n')
fprintf('error in labor mkt equilibrium  = %9.3e \n',  norm(sum(alpha, 2) - 1));
fprintf('error in wage markdowns         = %9.3e \n',  norm(err_delta(:)));
fprintf('\n')


%% Aggregation
% Aggregates in economy with markups and markdowns

MuDeltas= (sum(omega ./ (mu .* delta), 2)).^(-1);                   % sectoral markup x markdown

pratio  = (N .* mu .* delta .* alpha ./ MuDeltas).^(1/(1-p.eta));   % price ratio: p_i(s)/p(s)
Zs      = ( 1/N * sum(z.^(-1-1/p.etaL) .* pratio.^(p.eta*(-1-1/p.etaL)), 2) ...
          ).^(-p.etaL/(p.etaL + 1));                                % sectoral TFP

% -------------------------------------------------------------------------
% the formula for Ws is wrong
% wratio  = (N * alpha).^(1/(1 + p.etaL));                          % wage ratio: W_i(s)/W(s)                   
wratio  =  pratio .* z./Zs .* MuDeltas ./ (mu.*delta);

Ws      = (MuDeltas ./ Zs) ./ ...
          (1/N * sum((mu.*delta.*wratio./z).^(1-p.eta), 2)).^(1/(1-p.eta)); 
                                                                    % Ws obtained by equating two equations for omega(s)
% -------------------------------------------------------------------------
Omegas  = 1/S * (MuDeltas .* Ws ./ Zs).^(1-p.theta);                % sectoral sales share
MuDelta = 1/sum(Omegas ./ MuDeltas);                                % Aggregate markup x mardown

Ps      = (Omegas*S).^(1/(1-p.theta));                              % Sectoral price
Z       = (1/S * sum((Omegas * S).^(-p.theta/(1-p.theta)) ./ Zs))^(-1); % Aggregate TFP

W       = Z / MuDelta;
L       = (W * Z^(-p.sigma))^(1/(p.sigma + p.phi));                 % INTRA-temporal tradeoff
Y       = Z * L;
C       = Y;
Pi      = Y - W*L;                                                  % Aggregate profits


% Efficient Allocations: compute corresponding objects above in the absence
% of markups

Zse     = ( sum(z.^(p.eta - 1), 2) ).^(1/(p.eta - 1));

Ze      = ( 1/S * sum(Zse.^(p.theta - 1)) ).^(1/(p.theta - 1));
We      = Ze;
Le      = We^(1/(p.sigma + p.phi)) * Ze^(-p.sigma/(p.sigma + p.phi));
Ye      = Ze * Le; 
Ce      = Ye;
Ve      = Ce^(1-p.sigma)/(1-p.sigma) - Le^(1+p.phi)/(1+p.phi);              % efficient welfare

tau     = ( (Ve + L^(1+p.phi)/(1+p.phi)) * (1-p.sigma)/(C^(1-p.sigma)) ...
          ) ^(1/(1-p.sigma))  - 1 ;                                         % consumption equivalent of the welfare loss

% Print aggregate effects of markup distortions
fprintf('\n');
fprintf('Benchmark (left), Efficient (middle), Log Diff (right) \n');

fprintf('\n');
fprintf('Agg Markup*Markdown                    = %9.3f            \n',  MuDelta);
fprintf('Agg Productivity                       = %9.3f %9.3f %9.3f\n',  [Z, Ze, log(Ze/Z)]);
fprintf('Consumption                            = %9.3f %9.3f %9.3f\n',  [C, Ce, log(Ce/C)]);
fprintf('Employment                             = %9.3f %9.3f %9.3f\n',  [L, Le, log(Le/L)]);
fprintf('Welfare loss                           = %9.3f            \n',  tau);
fprintf('\n');

% FOLLOWING CODE HAS NOT BEEN MODIFIED FROM PS11 --------------------------
% Compute cost-weighted distribution of markups (via percentiles)

weight   = (z ./ Zs).^(p.eta - 1) .* (mu./Mus).^(-p.eta);                % labor share of a given firm

fprintf('\n')
fprintf('10 p.c Markup (cost weighted)          = %9.3f \n',  wprctile(mu(:), 10, 1/S*weight(:)));
fprintf('25 p.c Markup (cost weighted)          = %9.3f \n',  wprctile(mu(:), 25, 1/S*weight(:)));
fprintf('50 p.c Markup (cost weighted)          = %9.3f \n',  wprctile(mu(:), 50, 1/S*weight(:)));
fprintf('75 p.c Markup (cost weighted)          = %9.3f \n',  wprctile(mu(:), 75, 1/S*weight(:)));
fprintf('90 p.c Markup (cost weighted)          = %9.3f \n',  wprctile(mu(:), 90, 1/S*weight(:)));
fprintf('95 p.c Markup (cost weighted)          = %9.3f \n',  wprctile(mu(:), 95, 1/S*weight(:)));
fprintf('\n')

% Compute Concentration Statistics
omega_sorted = sort(omega, 1, 'descend');                                               % sort omega along each column (i.e. within sectors)
num_rows     = size(omega, 1);
row4         = min(4, num_rows); 

top1         = omega_sorted(1, :);                                                       % max mkt share per row (i.e. within any sector)
top4         = sum(omega_sorted(1:row4, :));
invHHI       = 1 ./ sum(omega.^2, 1);

fprintf('\n')
fprintf('median largest sales share             = %9.3f \n',  median(top1  ));
fprintf('median top 4 sales share               = %9.3f \n',  median(top4  ));
fprintf('inverse herfhindhal (goods)            = %9.3f \n',  median(invHHI));    
fprintf('\n')
% fprintf('median largest wage bill share         = %9.3f \n',  median(...  ));
% fprintf('median top 4 wage bill share           = %9.3f \n',  median(...  ));
% fprintf('inverse herfhindhal (labor)            = %9.3f \n',  median(...));    
% fprintf('\n')
% -------------------------------------------------------------------------



