import "trmodel"

module equilibrium_module (R:real) (trm:tr_model with t = R.t) : equilibrium_module with t = R.t {
    
    type t = R.t
    type mp [n][c][Ax][ns][nd] = trm.mp [n][c][Ax][ns][nd]
    type transition[ns] = trm.transition [ns]

    type sol [n][c][Ax][ns][nd][np] = 
        {p             : [np]t,         -- prices
         ed            : [np]t,         -- excess demand vector
         F             : transition[ns] -- struct with physical transition matrices
         ev_tau        : [n][ns]t       -- expected values for all types
         ccp_tau       : [n][ns][nd]t   -- conditional choice probabilities for all choices
         ccp_scrap_tau : [n][ns]t       -- scrapping probabilities dependent on not having a car
         ctp_tau       : [n][ns][ns]t   -- total transition matrix delta_tau*F.notrade
         delta_tau     : [n][ns][ns]t   -- trade transition probability matrices for all consumer types
         h_tau         : [n][ns]t       -- post-trade equilibrium holdings for all consumer types
        }
}