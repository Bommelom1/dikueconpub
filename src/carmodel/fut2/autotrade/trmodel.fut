
import "../lib/github.com/diku-dk/sparse/mono"

-- ********************
-- Car model parameters
-- ********************

module type trmodel = {

  -- basic real type
  type t

  val nanmax [n] : [n]t -> t
  val nansum [n] : [n]t -> t
  val logsum [n] : (sigma:t) -> [n]t -> t

  val nanmax2 : t -> t -> t
  val nansum2 : t -> t -> t
  val logsum2 : (sigma:t) -> t -> t -> t

  -- model parameters
  type mp [n][c][Ax][ns][nd]
  val mk : (n:i64) -> (c:i64) -> (Ax:i64) -> ?[ns][nd].mp[n][c][Ax][ns][nd]

  val set_newprices  [n][c][Ax][ns][nd] : mp[n][c][Ax][ns][nd] -> [c]t ->
                                          mp[n][c][Ax][ns][nd]

  val set_pscrap     [n][c][Ax][ns][nd] : mp[n][c][Ax][ns][nd] -> [c]t ->
                                          mp[n][c][Ax][ns][nd]

  val set_u_0        [n][c][Ax][ns][nd] : [n][c]t -> mp[n][c][Ax][ns][nd] ->
                                          mp[n][c][Ax][ns][nd]

  type consumertype = i64

  -- car prices (for new and old cars)
  type prices [c][Ax] = [Ax][c]t   -- c: ncartypes, Ax: maxage

  val simple_prices  [n][c][Ax][ns][nd] : mp[n][c][Ax][ns][nd] -> t -> prices[c][Ax]

  -- utilities of transitions
  type utility [ns][nd] = [ns][nd]t
  val utility        [n][c][Ax][ns][nd] : mp[n][c][Ax][ns][nd] -> prices[c][Ax]
                                          -> consumertype -> utility[ns][nd]

  type~ transition [ns] -- = {trade:sp.mat[ns][ns], notrade:sp.mat[ns][ns]}
  val age_transition [n][c][Ax][ns][nd] : mp[n][c][Ax][ns][nd] -> transition[ns]

  -- expected values
  type ev[ns] = [ns]t
  val ev0            [n][c][Ax][ns][nd] : mp [n][c][Ax][ns][nd] -> ev[ns]
  val vec_ev         [n][c][Ax][ns][nd] : mp [n][c][Ax][ns][nd] -> ev[ns] -> [ns]t

  -- the bellman equation
  val bellman        [n][c][Ax][ns][nd] : mp[n][c][Ax][ns][nd] -> utility[ns][nd]
                                          -> transition[ns] -> ev[ns] -> ev[ns]
  val bellmanJ       [n][c][Ax][ns][nd] : mp[n][c][Ax][ns][nd] -> utility[ns][nd]
                                          -> transition[ns] -> ev[ns]
					  -> (ev[ns], [ns][ns]t)

  -- for testing:
  val transition_notrade [ns] : transition[ns] -> [ns][ns]t
  val transition_trade   [ns] : transition[ns] -> [ns][ns]t

  type car = {cartype:i64, age:i64}
  type state = #NoCar | #Car car
  type decision = #Keep | #Purge | #Trade car

  val carprice_sell [n][c][Ax][ns][nd] : mp[n][c][Ax][ns][nd] -> prices[c][Ax] -> state -> t
  val carprice_buy [n][c][Ax][ns][nd] : mp[n][c][Ax][ns][nd] -> prices[c][Ax] -> decision -> t

  val ev_scrap [n][c][Ax][ns][nd] : mp[n][c][Ax][ns][nd] -> prices[c][Ax]
                                    -> consumertype -> state -> t

}

module trmodel (R:real) : trmodel with t = R.t = {

  module sparse = mk_mono R  -- or mk_compressed (mono sparse row representation is more efficient)
  module sp = sparse.sr

  type t = R.t

  type mp [n][c][Ax][ns][nd] =       -- n: ntypes, c: ncartypes, Ax: maxage
     {bet             : t,         -- discount factor
      mum             : [n]t,      -- marginal utility of money
      tw              : [n]t,      -- consumer type distribution; must sum to 1
      pnew            : [c]t,      -- car specific abar - we should tighten this
      pscrap          : [c]t,
      psych_transcost : [n]t,
      u_0             : [n][c]t,
      u_a             : [n][c]t,
      sigma           : t,
      sigma_s         : t,
      ns              : [ns](),
      nd              : [nd](),
      maxage          : [Ax]()
     }

  def mk (n:i64) (c:i64) (Ax:i64) : ?[ns][nd].mp[n][c][Ax][ns][nd] =
    {bet             = R.f32 0.95,
     mum             = replicate n (R.f32 0.1),
     tw              = replicate n R.(i64 1 / i64 n),
     pnew            = replicate c (R.i32 200),  -- kr 200.000
     pscrap          = replicate c (R.f32 1),    -- kr       1
     psych_transcost = replicate n (R.i32 0),
     u_0             = replicate n (replicate c (R.f32 6)),
     u_a             = replicate n (replicate c (R.f32(-0.5))),
     sigma           = R.f32 1,
     sigma_s         = R.f32 0.5,
     ns              = replicate (c*Ax+1) (),
     nd              = replicate (c*Ax+2) (),
     maxage          = replicate Ax ()
     }

  -- some utilities
  def mapi 'a 'b [n] (f: i64 -> a -> b) (xs:[n]a) : [n]b =
    map2 f (iota n) xs

  -- some setters
  def set_newprices [n][c][Ax][ns][nd]
                    (mp:mp [n][c][Ax][ns][nd])
                    (pnew:[c]t) : mp [n][c][Ax][ns][nd] =
    mp with pnew=pnew

  def set_pscrap [n][c][Ax][ns][nd]
                 (mp:mp [n][c][Ax][ns][nd])
                 (pscrap:[c]t) : mp [n][c][Ax][ns][nd] =
    mp with pscrap=pscrap

  def set_u_0 [n][c][Ax][ns][nd] (u_0:[n][c]t)
              (mp:mp [n][c][Ax][ns][nd]) : mp [n][c][Ax][ns][nd] =
    mp with u_0=u_0

  type consumertype = i64

  type prices [c][Ax] = [Ax][c]t   -- c: ncartypes

  def simple_prices [n][c][Ax][ns][nd] (mp:mp [n][c][Ax][ns][nd]) (r:t) : prices [c][Ax] =
    let is = ([R.i32 1] ++ scan (R.*) (R.i32 1) (replicate (Ax-1) r)) :> [Ax]t
    in transpose (map (\p -> map (R.* p) is) mp.pnew)

  type car = {cartype: i64, age: i64}
  type state = #NoCar | #Car car
  type decision = #Keep | #Purge | #Trade car

  def carprice_sell [n][c][Ax][ns][nd] (mp:mp[n][c][Ax][ns][nd]) (p:prices[c][Ax]) (s:state) : t =
    match s
    case #NoCar -> R.f32 0
    case #Car {cartype,age} ->
      assert (age >= 0 && age <= Ax)
      (if age == Ax then mp.pscrap[cartype]
       else p[age,cartype])

  def carprice_buy [n][c][Ax][ns][nd] (_mp:mp[n][c][Ax][ns][nd]) (p:prices[c][Ax]) (d:decision) : t =
    match d
    case #Keep -> R.nan
    case #Purge -> R.f32 0
    case #Trade {cartype,age} ->
      assert (age >= 0 && age < Ax)  -- notice age < Ax!!
       (p[age,cartype])

  def u_car [n][c][Ax][ns][nd] (mp:mp[n][c][Ax][ns][nd]) (tau:consumertype) ({cartype,age}: car) : t =
    let () = assert (tau >= 0 && tau < n) ()
    in if age == Ax then R.nan
       else R.(mp.u_0[tau,cartype] + mp.u_a[tau,cartype] * i64 age)

  def nanmax2 (x:t) (y:t) : t =
    R.(if isnan x
       then y else if isnan y
              then x else max x y)

  def nanmax (xs:[]t) : t =
    map (\x -> if R.isnan x then R.lowest else x) xs |>
    reduce R.max R.lowest

  def nansum2 (x:t) (y:t) : t =
    R.(if isnan x
       then y else if isnan y
              then x else x+y)

  def nansum (xs:[]t) : t =
    map (\x -> if R.isnan x then R.i32 0 else x) xs |>
    reduce (R.+) (R.i32 0)

  def logsum [n] (sigma:t) (v:[n]t) : t =
    let maxv = nanmax v
    in if R.(sigma == i32 0) then maxv
       else R.(map ((\x -> x - maxv) >-> (/sigma) >-> exp) v |>
               nansum |> log |> (*sigma) |> (+ maxv))

  def logsum2 (sigma:t) (x:t) (y:t) : t =
    let maxv = nanmax2 x y
    in R.(if sigma == i32 0 then maxv
          else let f v = v - maxv |> (/sigma) |> exp
               in nansum2 (f x) (f y) |> log |> (*sigma) |> (+ maxv)
         )

  def ev_scrap [n][c][Ax][ns][nd] (mp:mp[n][c][Ax][ns][nd]) (p:prices[c][Ax])
                                  (tau:consumertype) (s:state) : t =
    let () = assert (tau >= 0 && tau < n) ()
    let default () =
      let pscrap =
        match s
        case #NoCar -> R.nan
        case #Car {cartype,age=_} -> mp.pscrap[cartype]
      let psell = carprice_sell mp p s
      in logsum2 mp.sigma_s (R.(mp.mum[tau] * (pscrap-psell))) (R.i32 0)
    in match s
       case #NoCar -> default()
       case #Car{age,cartype=_} -> if age == Ax then R.i32 0
                                   else default()
  ------------------------
  -- Definition of utility
  ------------------------

  def util [n][c][Ax][ns][nd] (mp:mp[n][c][Ax][ns][nd]) (p:prices[c][Ax])
                              (tau:consumertype) (s:state) (d:decision) : t =
    assert (tau >= 0 && tau < n)
    (match d
     case #Keep ->
       (match s
        case #NoCar -> R.nan
        case #Car car -> u_car mp tau car
       )
     case #Purge ->
       (match s
        case #NoCar -> R.i32 0
        case #Car _ ->
          R.(mp.mum[tau] * carprice_sell mp p s +
             ev_scrap mp p tau s)
       )
     case #Trade newcar ->
       (match s
        case #NoCar ->
          R.(u_car mp tau newcar -
             mp.mum[tau] * carprice_buy mp p d -
             mp.psych_transcost[tau])
        case #Car _ ->
          R.(u_car mp tau newcar -
             mp.mum[tau] * carprice_buy mp p d -
             mp.psych_transcost[tau] +
             mp.mum[tau] * carprice_sell mp p s +
             ev_scrap mp p tau s)
       )
    )

  type utility [ns][nd] = [ns][nd]t

  def utility [n][c][Ax][ns][nd] (mp:mp[n][c][Ax][ns][nd]) (p:prices[c][Ax])
                                 (tau:consumertype) : utility[ns][nd] =
    -- notice: in decisions, car ages range from 0 to Ax-1 whereas
    -- in states, car ages range from 1 to Ax...
    let utils_keep : [ns]t =
      (flatten (
         tabulate_2d c Ax
           (\j a -> util mp p tau (#Car{cartype=j,age=a+1}) #Keep)
       ) ++ [util mp p tau #NoCar #Keep]
      ) :> [ns]t
    let utils_car : [][ns]t =
      flatten (
        tabulate_2d c Ax
          (\j' a' ->
             (flatten
               (tabulate_2d c Ax
                  (\j a ->
                    util mp p tau (#Car{cartype=j,age=a+1}) (#Trade{cartype=j',age=a'})))
              ++ [util mp p tau (#NoCar) (#Trade{cartype=j',age=a'})]) :> [ns]t)
      )
    let utils_purge : [ns]t =
      (flatten (
         tabulate_2d c Ax
           (\j a -> util mp p tau (#Car{cartype=j,age=a+1}) #Purge)
       ) ++ [util mp p tau #NoCar #Purge]
      ) :> [ns]t
    in transpose (([utils_keep] ++ utils_car ++ [utils_purge]) :> [nd][ns]t)

  type~ transition [ns] = {trade:sp.mat[ns][ns], notrade:sp.mat[ns][ns]}

  def age_transition [n][c][Ax][ns][nd] (_:mp[n][c][Ax][ns][nd]) : transition[ns] =
    let next_keep : [ns]i64 =
      (tabulate_2d c Ax (\j a ->
                           if a == Ax-1
                           then j*Ax -- old car (clunker) replaced with new car of same car type
                           else j*Ax+a+1)
       |> flatten
       |> (++[ns-1]) -- the keep state
      ) :> [ns]i64
    in {trade = sp.eye ns ns,
        notrade = sp.sparse ns ns (zip3 (iota ns) next_keep (replicate ns (R.i32 1)))}

  type ev[ns] = [ns]t
  def ev0 [n][c][Ax][ns][nd] (_:mp [n][c][Ax][ns][nd]) : ev[ns] =
    replicate ns (R.i32 0)

  def vec_ev [n][c][Ax][ns][nd] (_:mp [n][c][Ax][ns][nd]) (v:ev[ns]) : [ns]t =
    v

  def bellman0 [n][c][Ax][ns][nd]
               (mp:mp[n][c][Ax][ns][nd])
               (utils: utility[ns][nd])
               (F:transition[ns])
              : ev[ns] -> (ev[ns], [ns][nd]t) =
    \(ev:ev[ns]) ->
      let keep : [ns]t =                     -- keep: column 0 (of nd)
        sp.smvm F.notrade ev |> map (R.* mp.bet)
        |> map2 (R.+) (utils[:,0])
        |> mapi (\i v -> if i == ns-1 then R.nan else v) -- no car

      let trade_idx_n = nd-2
      let trade : [ns][trade_idx_n]t =                  -- trade: column 1 to nd-2 (of nd)
        let w : [ns]t =
          ev                                            -- was: sp.smvm F.trade ev
          |> map (R.* mp.bet)
        let w' : [ns][trade_idx_n]t =
          w[0:trade_idx_n] |> replicate ns
        let utils_trade : [ns][trade_idx_n]t =
          (utils[:,1:nd-1] :> [ns][trade_idx_n]t)
        in map2 (map2 (R.+)) w' utils_trade

      let purge : [ns]t =                   -- purge: column nd-1 (of nd)
        let s = mp.bet R.* ev[ns-1]
        in map (R.+ s) (utils[:,nd-1])

      -- let v : [ns][nd]t =
      --   ([keep] ++ trade ++ [purge] :> [nd][ns]t)
      --   |> transpose

      let v : *[ns][nd]t = tabulate_2d ns nd (\ _ _ -> R.i64 0)

      let v[:,0] = keep
      let v[:,1:trade_idx_n+1] = trade
      let v[:,nd-1] = purge

--      let v =
--        map2 (map2 (\u v -> if R.isnan u then u else v)) utils v

      in (map (logsum mp.sigma) v, v)

  def bellman mp utils F ev =
    bellman0 mp utils F ev |> (.0)

  def trade_transition [n][c][Ax][ns][nd]
                       (_mp:mp[n][c][Ax][ns][nd])
                       (ccp: [ns][nd]t) : [ns][ns]t =
    let deltaKdiag : [ns]t = map (\x -> x[0]) ccp -- keep: column 0
    let deltaT = map (\(x:[nd]t) ->
		        let vs : [c][Ax]t = unflatten (x[1:ns] :> [c*Ax]t)
		        let vs = map (rotate 1) vs |> flatten
			in (vs ++ [x[nd-1]]) :> [ns]t
		     ) ccp
    in mapi (\i -> mapi (\j e -> if i == j then e R.+ deltaKdiag[i]
				 else e)
	    ) deltaT

  def bellmanJ [n][c][Ax][ns][nd]
               (mp:mp[n][c][Ax][ns][nd])
               (utils: utility[ns][nd])
	       (F:transition[ns])
	       (ev:ev[ns]) : (ev[ns], [ns][ns]t) =
    let (ev1:ev[ns], v:[ns][nd]t) = bellman0 mp utils F ev
    let ccp : [ns][nd]t =
      map2 (\r x -> map (R.exp <-< (R./ mp.sigma) <-< (R.- x)) r) v ev1
    let ccp = map (map (\x -> if R.isnan x then R.i64 0 else x)) ccp
    let delta : [ns][ns]t = trade_transition mp ccp
    --let dev = sp.dmsmm delta F.notrade |> map (map (R.* mp.bet))
    let dev = sp.scale mp.bet F.notrade |> sp.dmsmm delta
    in (ev1,dev)

  def transition_notrade [ns] (tr:transition[ns]) : [ns][ns]t =   -- for testing
    sp.dense (tr.notrade)

  def transition_trade [ns] (tr:transition[ns]) : [ns][ns]t =     -- for testing
    sp.dense (tr.trade)
}
