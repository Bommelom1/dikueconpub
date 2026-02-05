%% Estimation with new data and endogenous scrappage 
% 2021-9-11

clc; 
close all; 
clear all; 
addpath('matlabinclude');
addpath('autotrade');

years = 1996:2008;   % sample period: max 1996:2009

%% -------------- Read data -------------------
% Read in data from .xlsx files or prev. stored .mat files
% Type "help data.setup" for a list of available datasets and instructions on how to add new data sets
datadir='data/8x4/'; 

% set up parameters and indexes consistent with loaded data
LOADMATFILES = true; % set to false first time you run this program to generate matfiles 
[mp, s, dta, car, demo, fuelp, est]=data.setup(datadir, years, LOADMATFILES); 

%  Estimate model outcomes (ccps, holding distributions ect. ) based on tabulated choice/state/type data
[sol_d] = estim.outcomes(mp, s, dta); 

save('results/estimation/estimation_data.mat', 'mp', 's', 'dta', 'car', 'demo', 'fuelp', 'est', 'sol_d')


%% STEP 1: Estimate restricted model

mp = trmodel.setparams(mp); 

mp.pnames = { 'u_0', 'u_a', 'u_even','mum', 'psych_transcost', 'psych_transcost_nocar', 'acc_0', 'acc_a', 'sigma_s'}; 
mp.pnames = {mp.pnames{:}, 'tc_sale', 'tc_sale_even'};

mp.p0=getfields(mp, mp.pnames); % set starting values to parameters in mp

mp.p0.u_0     = repcell({3.5}, mp.ntypes, mp.ncartypes);
mp.p0.u_a     = repcell({-0.5}, mp.ntypes, mp.ncartypes);
mp.p0.u_a_sq  = repcell({0}, 1, mp.ncartypes);
mp.p0.u_even  = repcell({0}, mp.ntypes, mp.ncartypes);
mp.p0.acc_0   = repcell({-5}, 1, mp.ncartypes); % logit transformed parameter 
mp.p0.acc_a   = repcell({0}, 1, mp.ncartypes); % logit transformed parameter 
mp.p0.mum     = repcell({0.1}, mp.ntypes,1);
mp.p0.psych_transcost = repcell({5}, mp.ntypes, 1); 
mp.p0.tc_sale=0;
mp.p0.sigma_s =1;

% fixed parameters 
mp.pscrap_notax=num2cell([mp.pnew_notax{:}]*.87^mp.abar_j0{1});  
mp.pscrap_notax={10};
mp.pscrap_notax=num2cell([mp.pnew_notax{:}].*car.beta(car.age==0)'.^ [mp.abar_j0{:}]);  
mp.pscrap_notax={-6};
mp.pscrap_notax=num2cell([mp.pnew_notax{:}].*car.beta(car.age==0)'.^ [mp.abar_j0{:}]);  

% update parameters
mp = trmodel.setparams(mp);

% estimate model
[mp_mle, theta_mle, Avar_mle, sol_mle, pvec0_mle] = estim.mle(mp, s, dta, {'bhhh','bhhh','bhhh'}, [20,20,10]);

save('results/estimation/mle_step1.mat', 'mp_mle', 'sol_mle', 'Avar_mle')

%% STEP 2: Run again
load('results/estimation/mle_step1.mat')
mp = mp_mle; 
mp.p0=getfields(mp, mp.pnames); % set starting values to parameters in mp
[mp_mle, theta_mle, Avar_mle, sol_mle, pvec0_mle] = estim.mle(mp, s, dta, {'bhhh','bhhh','bhhh'}, [10,10,10]);

save('results/estimation/mle_step2.mat', 'mp_mle', 'sol_mle', 'Avar_mle')



