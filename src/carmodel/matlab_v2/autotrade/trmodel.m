% The main static class for the AUTOMOBILE TRADE model.
% Contains all needed model parts.
% 
% November 2021

classdef trmodel
methods (Static)
  function mp = default_params()
    % Default parameters used for examples
    % SYNTAX :
    %   [mp] =setparams.default_params
    %
    % See also trmodel.setparams to get list of activate parameters used in trmodel (all set to)
 

    mp.abar = 25;                 % car specific value of abar 
    mp.tw=[.5 .5]';                 % distribution of the types in the population (tw vector must sum to 1)
 
    mp.sigma=1.;          % scale value of extreme value taste shocks in consumer problem at top level of the choice tree
    mp.sigma_s=0.5; 

    mp.mum = .1;            % marginal utility of money
    mp.psych_transcost = 0; 
    mp.u_0=6;               % intercept in utility function (ntypes x ncartypes cell) 
    mp.u_a=-.5;             % slope coefficient on car-age in utility function (ntypes x ncartypes cell) 
        
    % Prices before taxes (prices after taxes are computed when solving the model)
    mp.pnew=200;            % new car prices (cartype specific)
    mp.pscrap=1;            % scrap  prices (cartype specific)
    mp.p_fuel=10.504;       % pricer of fuel (DKK/liter) from average fuel price in 2008 from our dataset    
  end % end of default_params

  function [mp] =setparams(mp0)
    % Standard parameters of the model
    % SYNTAX :
    %   [mp] = trmodel.setparams       set all parameters to zeros
    %   [mp] = trmodel.setparams(mp0)  update patameters with parameters in mp0
    %
    % NOTE: When adding new parameters to trmodel, they must be added here 
    % DO NOT ENTER VALUES OF PARAMETERS - JUST ZEROS (or values that turn off the parameter)
    % See also setparams.default to get a default parameters

    % ******************************************
    % consumer types and car types
    % ******************************************
    % Some parameters below are consumer type specific (in rows) and car type specific (in columns) 
    % If number of entries is smaller than ntypes, last entry is repeated)

    mp.abar = 25;  % car specific value of abar 
    mp.tw=1;          % distribution of the types in the population (tw vector must sum to 1)
    mp.ncartypes=1;   % number of car types 
    
    % ******************************************
    % discount factor (beta)
    % ******************************************
    mp.bet=.95;  

    % ******************************************
    % parameters of GEV distribution (theta)
    % ******************************************
    mp.sigma=1;        % scale value of extreme value taste shocks in consumer problem at top level of the choice tree
    mp.sigma_s=1;       % the extreme value scale parameter for the idiosyncratic extreme value shocks affecting the sale vs scrappage decision
 
    % ******************************************
    % utility patameters (theta)
    % ******************************************
    mp.mum = 1;       % marginal utility of money
    
    % transactions costs parameters 
    mp.psych_transcost      = 0;      % utility cost of bying a car

    % car utility parameters (see also trmodel.u_car) 
    % Reduced form car utility mp.u_0{tp, car}+mp.u_a{tp, car}*car_age + mp.u_a_sq{tp, car}*car_age.^2
    mp.u_0=0;           % intercept in utility function (ntypes x ncartypes cell) 
    mp.u_a=0;           % slope coefficient on car-age in utility function (ntypes x ncartypes cell) 

    % ******************************************
    % Prices 
    % if car or consumer type parameters have size smaller than ncartypes and ntypes, the last element is repeated.
    % ******************************************
    mp.pnew=100;    % new car prices (cartype specific)
    mp.pscrap=0;    % scrapcar prices (cartype specific)
    mp.p_fuel=5;    % pricer of fuel (DKK/liter) from average fuel price before tax
    
    % ********************************************************************************
    % update parameters
    % ********************************************************************************
    mp.ntypes=numel(mp.tw);  % number of types of consumers

    % update mp with mp0 (default values are over written by input values mp0)
    if nargin>0
      mp=combinestructs(mp, mp0);
    end

    % update endogenous parameters
    mp=trmodel.update_mp(mp);
  end % end of setparams

  function mp = update_mp(mp)
    % update_mp: Updates mp structs so that dimensions of fields in 
    % 
    % Purpose: If number of elements in fields is smaller than the required dimension, 
    % then last element is repeated using the function extend 
    % If the required dimension is smaller that the number of elements in a field, 
    % only the first required elements are returned so that mp are consistent with model. 

    % Parameters that vary over car types, j
    param_j       = {'pnew', 'pscrap'};
    for k=1:numel(param_j)
      mp.(param_j{k}) =extend(mp.(param_j{k}),1, mp.ncartypes);
    end  


    % Parameters that vary over consumer types, tau
    param_tau     = {'mum','psych_transcost'};
    for k=1:numel(param_tau)
      mp.(param_tau{k}) =extend(mp.(param_tau{k}),mp.ntypes,1);
    end  

    % Parameters that vary over consumer types, tau and car types, j
    param_tau_j   = {'u_0','u_a'};
    for k=1:numel(param_tau_j)
      mp.(param_tau_j{k}) =extend(mp.(param_tau_j{k}),mp.ntypes, mp.ncartypes);
    end  
  
    if (mp.ntypes == 1)
       mp.tw=1;
    end

    if (abs(sum(mp.tw)-1.0) > 1e-12)
       fprintf('Error in trmodel.setup: mp.tw vector does not sum to 1\n');
    end   
  end % end of update_mp

  function [s] = index(mp)
    % This function set up the state space for the trmodel. 
    % 
    % SYNTAX: [s] = trmodel.index(mp)
    % 
    % INPUT: 
    %   mp:     structure with model parameters (see trmodel.setparams)
    %
    % OUTPUT:  
    %   s:  a struct with the following fields:
    %
    %     ns: number of states
    %     nd: number of decisions
    %     id: struct that holds decision indexes groups of decisions. 
    %     is: struct that holds state indexes for groups of states
    %     tr: struct with sub-indexes for transitions
    %
    %  Example: When mp.ncartypes=2 and mp.abar=5 
    %  s=trmodel.index(mp) gives
    % s = 
    %
    %   struct with fields:
    %
    %     ns: 11
    %     nd: 12
    %     is: [1x1 struct]
    %     id: [1x1 struct]
    %     tr: [1x1 struct]
    % 
    %  where
    % 
    %   s.id has the following fields
    %                keep: 1
    %               trade: {[2 3 ... 17]  [18 19 ... 39]}
    %          trade_used: {[3 4 ...  17]  [19 20 ...  39]}
    %           trade_new: {[2]  [18]}
    %                 age: [NaN 0 1 ...  15 0 1 2 ...  21 NaN]
    %               purge: 40
    %    
    %   s.is has the following fields
    %      keep: 1
    %     trade: [10x2 double]
    % trade_vec: [2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21]
    %     purge: 22
    %       age: [NaN 0 1 2 3 4 5 6 7 8 9 0 1 2 3 4 5 6 7 8 9 NaN]
    %  
    % 
    %  usage: 
    %   s.id are choice indexes that for example used for indexing columns in matrices 
    %   such as ccp_tau, and decision specific utility
    %
    %   s.is are state indexes that for example used for indexing rows in matrices 
    %   such as ccp_tau, ev_tau, q, q_tau{tau}.
    %
    %   For example, ccp_tau{tau}(s.is.car(1:end-1,j), s.id.keep) gives the vector of probabilities 
    %   of keeping a car of type j with ages 1,2,..,abar-1. The age of such a car can also be accessed by 
    %   s.is.age(s.is.car(1:end-1,j))  
    %
    % Some example descriptions of sub-indices: 
    % 
    % s.is.car (mp.abar x mp.ncartypes array): 
    %   s.is.car(:,j) is a matrix holding the indexes in the overall ordering of car states for car type j - except the no car state. 
    % s.is.nocar (scalar): 
    %   index for the state of having no car
    % s.is.age (s.ns vector): 
    %   s.is.age vector of car ages corresponding to state indexes. Age of no car is coded as NaN
  
    s.ns=mp.ncartypes*mp.abar+1;   % number of states      (add one for no car state)
    s.nd=mp.ncartypes*mp.abar+2;   % number of decisions   (add one for keep and purge)
    
    % state indexes
    s.is.car_vec  = 1:mp.ncartypes*mp.abar;
    s.is.car=reshape(s.is.car_vec,mp.abar, mp.ncartypes); 
    s.is.nocar=s.ns;  % nocar state
    s.is.age=[repmat((1:mp.abar), 1, mp.ncartypes) nan];       % not possible to enter state with new car

    % decision indexes
    s.id.keep=1;              % keep decision
    s.id.trade=s.is.car + 1;  % add 1 to account for keep being the first index
    s.id.trade_vec = reshape(s.id.trade, 1, numel(s.id.trade));
    s.id.purge=s.ns+1;        % purge decision
    s.id.age=[nan repmat((0:mp.abar-1), 1, mp.ncartypes) nan]; % not possible to buy clunker

    % index for trade transition matrix corresponding to columns in CCP matrix (see definition of deltaT in paper)
    s.tr.deltaT=[reshape([s.id.trade(2:end,:); s.id.trade(1,:)],1, s.ns-1) s.id.purge];

    % index for next period car state - when keeping (used in bellman) 
    s.tr.next_keep=(1:s.ns)+1;                      % state after keep (current car ages one period) 
    s.tr.next_keep(s.is.car(end,:))=s.is.car(1,:);  % clunkers replaced by new cars
    s.tr.next_keep(s.is.nocar)=s.is.nocar;          % no car - stays no car
  end % end of index

  function [u, ev_scrap, ccp_scrap]=utility(mp, s, tau, prices)
    %   utility:  Function to compute utility for each consumer in each state and for each decision 
    %             and its derivative wrt to used car prices. 
    %
    %   SYNTAX:     [u, ev_scrap, ccp_scrap]=trmodel.utility(mp, s, prices)
    %
    %   INPUTS:
    %     mp          structure with model parameters (see trmodel.setparams)
    %     s:          structure with indexes for dimensions. see trmodel.index
    %     tau:        household type 
    %     prices      prices is a 2 dimensional array with mp.ncartypes columns 
    %
    %
    %   OUTPUTS          
    %   u:      Matrix that store the utility
    %           received from trading for a car of type j and age a for a car of type j' and age d, or "purge"
    %           trade the current car for the outside good, or to keep the car, for each consumer type.
    %
    %           u has s.ns rows and s.nd columns and rows is indexed by s.is and columns by s.id. 
    %
    %           Rows correspond to the different states 
    %           i.e the type/age of the currently held car or the state of no car. 
    %           The last row corresponds to the state of having no outside good. ()
    %
    %           The columns correspond to the different possible decisions a consumer has.
    %           The different possible decisions are: keep, trade and purge. 
    %
    %           The first column is the utility of the decision to keep (s.id.keep=1) 
    %           for each possible age of car (or outside good in last column)
    %
    %           NOTE: The utility of keep contain a missing value in in rows where the age equals the scrappage age (is.clunker{j})
    %           for each car type, since a consumer is not allowed to keep a car once it reaches the scrappage age, and for a consumer
    %           who has no car, the law row of the u_tau provides the utility of continuing to stay with the outside good (keeping outside good)
    % 
    %           The last column, s.id.purge=s.nd, is the utility of the decision to purge 
    %            (i.e. choose the outside good, i.e. move to having no car)
    %
    %           The columns in between, from 2 to s.nd-1 or s.id.trade_vec, correspond to the decision to trade 
    %           for a car of age a, a=0,...,abar-1 for each of the ncartypes of cars
    % 
    %  ev_scrap  the extra option value from being able to scrap instead of sell a used car (only non-zero if mp.es=1)
    %            This is the s.ns x 1 vector of integrated value of scrap choice relative to selling choice
    %            i.e. logsum(mp.mum(tau)*(pscrap-psell), 0)  
    %  ccp_scrap: 
    %           s.ns x 1 vector of scrapping probabilities (conditional on not keeping)
    %
    %           Gillingham, Iskhakov, Munk-Nielsen, Rust and Schjerning, September 2020

    % net sales price after transactions cost 
    [psell]=trmodel.psell(mp, s, prices); 

    % net purchase price after transactions cost 
    [pbuy]=trmodel.pbuy(mp, s, prices); 
  
    u=nan(s.ns,s.nd); % utility: number of states by number of decisions (see state_decision_index for information about of indexing)

    % utility of scrap option (when selling)
    [ccp_scrap, ev_scrap]=trmodel.p_scrap(mp, s, prices, tau);

    for j=1:mp.ncartypes
      % utility of keeping (when feasible)
      u(s.is.car(1:end-1,j), s.id.keep) =  trmodel.u_car(mp, s.is.age(s.is.car(1:end-1,j)), tau, j);
      
      % utility of trading 
      u(: , s.id.trade(:,j))           =  trmodel.u_car(mp, s.id.age(s.id.trade(:,j)), tau, j) ...
                                      + ev_scrap-mp.mum(tau)*(pbuy(:,s.id.trade(:,j))-psell)-mp.psych_transcost(tau);
    end % end loop over car-types

    % utility of purging (i.e. selling/scrapping the current car and not replacing it, i.e. switching to the outside good)
    u(:, s.id.purge)                  =  mp.mum(tau)*psell + ev_scrap;      
    u(s.is.nocar , s.id.keep)         =  nan; % not possible to keep no car
  end % end of utility

  function [ccp_scrap, ev_scrap]=p_scrap(mp, s, prices, tau)
    %   p_scrap:  Function to compute expected utility of the option to scrap rather than sell  
    %             for each consumer in each car state, the probability of scrapping,  and the
    %             derivatives with respect to mum{t} and the sell-side transaction cost parameters
    %
    %   SYNTAX:   [ccp_scrap, ev_scrap]=trmodel.p_scrap(mp, s, psell, dpsell, tau)

    psell=trmodel.psell(mp, s, prices); 

    pscrap=nan(s.ns,1);
    for j=1:mp.ncartypes
      pscrap(s.is.car(1:end-1,j))=mp.pscrap(j);
    end

    ev_scrap = logsum([mp.mum(tau)*(pscrap-psell), zeros(s.ns,1)], mp.sigma_s); 
    ccp_scrap=1-exp(-(ev_scrap)/mp.sigma_s);

    for j=1:mp.ncartypes
      ccp_scrap(s.is.car(end,j))=1;
      ev_scrap(s.is.car(end,j))=0;
    end
  end % end of p_scrap
  
  function [uv]=u_car(mp, car_age, tau, j)
    % car utility
    % INPUTS: 
    %     j_age: vector 
    %     tau: household type, scalar index
    %     j: car type, scalar index 

    uv=mp.u_0(tau, j)+mp.u_a(tau, j)*car_age;
  end % end of u

  function [pbuy] = pbuy(mp, s, prices)         
    pbuy=nan(1, s.nd); 
    for j=1:mp.ncartypes
      pbuy(s.id.trade(1,j))=  mp.pnew(j);  
      pbuy(s.id.trade(2:end,j))= prices(:,j);  
      pbuy(s.id.purge)=0;
    end
  end

  function [psell] = psell(mp, s, prices)     
    % psell is the selling price net of seller-side transactions costs
    psell=nan(s.ns,1);
    for j=1:mp.ncartypes
      psell(s.is.car(1:end-1,j))=prices(:,j);
      psell(s.is.car(end,j))=mp.pscrap(j);
      psell(s.is.nocar)=0;

    end
  end % end of trade_cost

  function [F] = age_transition(mp, s)
    % age_transition.m: age transition probability matrices 
    % 
    %  SYNTAX:   [F]= trmodel.age_transition(mp, s)
    %
    % INPUTS: 
    %   mp:     structure with model parameters
    %   s:      structure with state and decision indexes (generated by trmodel.states) 
    %
    % OUTPUS: 
    %   F:      structure with age transition matrices 
    %
    %   The fields of F are block-diagonal matrices with a block for each car and 1 for no car  
    %     F.notrade:   s.ns x s.ns extended physical transition matrix (Similar to Q matrix in paper)
    %     F.trade:     s.ns x s.ns state transition matrix conditional on trading a car 

    F.notrade = sparse((1:s.ns),s.tr.next_keep,1,s.ns,s.ns);
    F.trade   = speye(s.ns);
  end

  function [delta] = trade_transition(mp, s, ccp, ccp_scrap)
    % Trade transition matrix
    deltaK=sparse(1:s.ns,1:s.ns,ccp(:,s.id.keep),s.ns,s.ns);

    % Trading transition probabilities
    deltaT     = ccp(:,s.tr.deltaT);

    % trade transition probability matrix
    delta=deltaT +  deltaK; 
  end % end of trade_transition

  function [ev1, ccp, dev] = bellman(mp, s, util, F, ev)
    % bellman.m:    Implements the Bellman operator ev1=Gamma(ev) for the consumer's problem of trading cars,
    %               to maximize discounted utility in the presence of a secondary market for cars
    %               This function returns the *expected value function* ev and is a fast vectorized
    %               program that also calculates choice probabilities implied by current expected value ev
    %               and the derivative of the Bellman operator
    %
    % INPUTS: 
    %   mp:     structure with model parameters
    %   s:      structure with state and decision indexes (generated by trmodel.states) 
    %   util:   utility matrix (n.s x n.s+1) (states in rows, decision in columns)
    %   Q:      s.ns x s.ns extended physical transition matrix (block-diagonal matrix with a block for each car and 1 for no car)
    %   F:      s.ns x s.ns state transition matrix conditional on trading a car (block-diagonal matrix with a block for each car and 1 for no car)
    %   ev:     s.ns x 1 vector of expected values 
    %
    % OUTPUS: 
    %   ev1:    s.ns  x 1 vector of updated expected values after evaluating the Bellman operator 
    %   ccp:    s.ns  x s.ns + 1 matrix of conditional choice probabilities  
    %   dev:    s.ns  x s.ns matrix of derivatives of Bellman operator
    % 
    %  Fixed point of Bellman can be found by calling dpsolver
    % 
    %  SYNTAX: 
    %     [ev,ccp,dev]= dpsolver.poly(@(ev) trmodel.bellman(mp, s, util, F, ev), ev0, mp, mp.bet);
    % 
    % See also:
    %   dpsolver.poly

    v=nan(s.ns,s.nd);        % choice specific value functions (states by decision)
    
    % Calculate the expected discounted utility of keeping car

    % Keeping is only valid for car ages  a=1,...abar-1 (i.e not in clunker state and nocar state)
    v(s.is.car_vec,s.id.keep) =  util(s.is.car_vec,s.id.keep) + mp.bet*F.notrade(s.is.car_vec,:)*ev; 

    % calculate the values for all options involving buying a car
    v(:, s.id.trade_vec)         =util(:,s.id.trade_vec) + (mp.bet*F.trade(s.is.car_vec,:)*ev)';   

    % calculate the values for the purge decision (ev of no car)
    v(:,s.id.purge)               =util(:,s.id.purge)   + mp.bet*ev(s.is.nocar);                     % calculate the values for the purge decision (or keep having
    v(isnan(util))=nan;
    ev1=logsum(v, mp.sigma); 

    if nargout>1
      ccp=exp((v-ev1)/mp.sigma);
      ccp(isnan(ccp))=0;  % restore the nans in the elements of ccp_cell corresponding to infeasible choices
    end

    if nargout>2
      %  compute trade transition probability matrix
      delta=trmodel.trade_transition(mp, s, ccp);

      % derivative of Bellman operator wrt ev
      dev=mp.bet*delta*F.notrade;
    end
  end % end of bellman

end %methods

end %class
