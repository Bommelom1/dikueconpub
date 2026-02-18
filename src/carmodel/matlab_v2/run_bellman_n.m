% This solves the bellman equation for a given set of prices and parameters
% Note that any calls require these addpaths 

% addpath('matlabinclude');
% addpath('autotrade');

% clear
% clc

% read default parameters
function [ev_sa, ev_poly]=run_bellman_n(n, c, abar, type, pri)

    
    if nargin<4
        type = "all";
    end
    if nargin<5
        pri = true;
    end

    mp0=trmodel.default_params;

    mp0.tw =  ones(n,1)/n;
    
    % Parameter adjustments (relative to defaults)
    mp0.ncartypes=c;
    mp0.abar=abar;
    mp0.ntypes=n;
    
    % populate mp with remaining fields and update model parameter dependencies
    mp = trmodel.setparams(mp0);
    
    % model indexing
    s = trmodel.index(mp);
    
    % initialize prices where to solve bellman
    p = mp.pnew.*0.85.^(1:mp.abar-1)'; 
    
    % pre-compute age transition matrices
    F = trmodel.age_transition(mp, s);
    
    ap=dpsolver.setup;
    if pri
        ap.printfxp=1;  % show output
    end
    ev0 = zeros(s.ns,1); % initial guess on expected value function

    ev_sa   = zeros(s.ns, mp.ntypes);
    ev_poly = zeros(s.ns, mp.ntypes);

    for t  = 1:mp.ntypes
       if pri
        fprintf('Solving model for sonsumer type, tau=%d\n ', t);
       end
        
        % pre-compute utility
        util = trmodel.utility(mp, s, t, p);   % update the current period, post-trade payoffs stored in the global util_tau
        
        % step 1: solve the DP problem for the price vector p for each consumer type
        bellman=@(ev) trmodel.bellman(mp, s, util, F, ev);
        
        if type~="poly"
            ev_sa(:,t)= dpsolver.sa(bellman, ev0, ap);
        end
        if type~="sa"
            ev_poly(:,t)= dpsolver.poly(bellman, ev0, ap, mp.bet);
        end
    end
end

% evaluate excess demand for a given price vector, p
% [ed, ded, sol]=equilibrium.edf(mp, s, p);

% solve for equilibrium price
% [sol]=equilibrium.solve(mp, s); % solve model in baseline

