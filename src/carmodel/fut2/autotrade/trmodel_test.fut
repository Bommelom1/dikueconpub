
import "trmodel"

module trm = trmodel f64

-- ==
-- entry: test_simple_prices
-- input { 2i64 [200f64,200f64] 3i64 }
-- output { [ [ 200f64, 200f64 ],
--            [ 170f64, 170f64 ],
--            [ 144.5f64, 144.5f64 ] ] }
-- input { 2i64 [170f64] 2i64 }
-- output { [ [ 170f64 ],
--            [ 144.5f64 ] ] }

entry test_simple_prices [c] (n:i64) (newprices:[c]f64) (Ax:i64) : trm.prices[c][Ax] =
  let [ns][nd] mp : trm.mp [n][c][Ax][ns][nd] = trm.mk n c Ax
  let mp = trm.set_newprices mp newprices
  in trm.simple_prices mp 0.85

-- ==
-- entry: test_notrade
-- input { 1i64 2i64 3i64 }
-- output { [ [ 0f64, 1f64, 0f64, 0f64, 0f64, 0f64, 0f64 ],
--            [ 0f64, 0f64, 1f64, 0f64, 0f64, 0f64, 0f64 ],
--            [ 1f64, 0f64, 0f64, 0f64, 0f64, 0f64, 0f64 ],
--            [ 0f64, 0f64, 0f64, 0f64, 1f64, 0f64, 0f64 ],
--            [ 0f64, 0f64, 0f64, 0f64, 0f64, 1f64, 0f64 ],
--            [ 0f64, 0f64, 0f64, 1f64, 0f64, 0f64, 0f64 ],
--            [ 0f64, 0f64, 0f64, 0f64, 0f64, 0f64, 1f64 ] ] }

entry test_notrade (n:i64) (c:i64) (Ax:i64) : [][]f64 =
  let [ns][nd] mp : trm.mp [n][c][Ax][ns][nd] = trm.mk n c Ax
  in trm.age_transition mp |> trm.transition_notrade

-- ==
-- entry: test_trade
-- input { 1i64 2i64 3i64 }
-- output { [ [ 1f64, 0f64, 0f64, 0f64, 0f64, 0f64, 0f64 ],
--            [ 0f64, 1f64, 0f64, 0f64, 0f64, 0f64, 0f64 ],
--            [ 0f64, 0f64, 1f64, 0f64, 0f64, 0f64, 0f64 ],
--            [ 0f64, 0f64, 0f64, 1f64, 0f64, 0f64, 0f64 ],
--            [ 0f64, 0f64, 0f64, 0f64, 1f64, 0f64, 0f64 ],
--            [ 0f64, 0f64, 0f64, 0f64, 0f64, 1f64, 0f64 ],
--            [ 0f64, 0f64, 0f64, 0f64, 0f64, 0f64, 1f64 ] ] }

entry test_trade (n:i64) (c:i64) (Ax:i64) : [][]f64 =
  let [ns][nd] mp : trm.mp [n][c][Ax][ns][nd] = trm.mk n c Ax
  in trm.age_transition mp |> trm.transition_trade

-- ==
-- entry: test_carprice_sell
-- input { 2i64 [100f64,100f64] 2i64 }
-- output { [85f64, 1f64, 85f64, 1f64, 0f64] }

entry test_carprice_sell [c] (n:i64) (newprices:[c]f64) (Ax:i64) : ?[ns].[ns]f64 =
  let [ns][nd] mp : trm.mp [n][c][Ax][ns][nd] = trm.mk n c Ax
  let mp = trm.set_newprices mp newprices
  let p = trm.simple_prices mp 0.85
  let ss = tabulate_2d c Ax (\j a -> #Car{cartype=j,age=a+1})
           |> flatten |> (++ [#NoCar])
  in map (trm.carprice_sell mp p) ss

-- ==
-- entry: test_ev_scrap
-- input { 2i64 [100f64,100f64] 2i64 }
-- output { [2.528264977846713e-8f64, 0f64, 2.528264977846713e-8f64, 0f64, 0f64] }

entry test_ev_scrap [c] (n:i64) (newprices:[c]f64) (Ax:i64) : ?[ns].[ns]f64 =
  let [ns][nd] mp : trm.mp [n][c][Ax][ns][nd] = trm.mk n c Ax
  let mp = trm.set_newprices mp newprices
  let p = trm.simple_prices mp 0.85
  let tau = 1
  let ss = tabulate_2d c Ax (\j a -> #Car {cartype=j,age=a+1})
           |> flatten |> (++ [#NoCar])
  in map (trm.ev_scrap mp p tau) ss

-- ==
-- entry: test_carprice_buy
-- input { 2i64 [100f64,100f64] 2i64 }
-- output { [f64.nan,100f64,85f64,100f64,85f64,0f64] }

entry test_carprice_buy [c] (n:i64) (newprices:[c]f64) (Ax:i64) : ?[nd].[nd]f64 =
  let [ns][nd] mp : trm.mp [n][c][Ax][ns][nd] = trm.mk n c Ax
  let mp = trm.set_newprices mp newprices
  let p = trm.simple_prices mp 0.85
  let ds = tabulate_2d c Ax (\j a -> #Trade{cartype=j,age=a})
           |> flatten |> (\tr -> [#Keep] ++ tr ++ [#Purge])
  in map (trm.carprice_buy mp p) ds


-- ==
-- entry: test_acc_prob
-- input { 2i64 2i64 3i64 }
-- output { [0.00004539786f64,0.00004539786f64,0.00004539786f64,0.00004539786f64,0.00004539786f64,0.00004539786f64,0f64] }

entry test_acc_prob (n:i64) (c:i64) (Ax:i64) : ?[ns].[ns]f64 =
  let [ns][nd] mp : trm.mp [n][c][Ax][ns][nd] = trm.mk n c Ax
  in trm.acc_prob mp

-- ==
-- entry: test_acc_prob_mat
-- input { 2i64 2i64 2i64 }
-- output { [ [0.0f64,0.00004539786f64,0f64,0f64,0f64],
--            [0.0f64,0.00004539786f64,0f64,0f64,0f64],
--            [0.0f64,0f64,0f64,0.00004539786f64,0f64],
--            [0.0f64,0f64,0f64,0.00004539786f64,0f64],
--            [0.0f64,0f64,0f64,0f64,0f64] ]}

entry test_acc_prob_mat (n:i64) (c:i64) (Ax:i64) : ?[ns].[ns][ns]f64 =
  let [ns][nd] mp : trm.mp [n][c][Ax][ns][nd] = trm.mk n c Ax
  let mat = trm.acc_prob_mat mp|> (.0)
  in trm.dense_acc_prob_mat mat

-- ==
-- entry: test_utility
-- input { 2i64 [100f64,100f64] 2i64 }
-- output { [ [5.5f64, 4.500000002930908f64, 5.50000002528265f64, 4.500000002930908f64,
--             5.50000002528265f64, 8.50000015194252f64],
--            [f64.nan, -3.900000147521496f64, -2.900000125169754f64, -3.900000147521496f64,
--             -2.900000125169754f64, 0.10000000149011612f64],
--            [5.5f64, 4.500000002930908f64, 5.50000002528265f64, 4.500000002930908f64,
--             5.50000002528265f64, 8.50000015194252f64],
--            [f64.nan, -3.900000147521496f64, -2.900000125169754f64, -3.900000147521496f64,
--             -2.900000125169754f64, 0.10000000149011612f64],
--            [f64.nan, -4.000000149011612f64, -3.00000012665987f64, -4.000000149011612f64,
--             -3.00000012665987f64, 0.0f64]]
-- }

entry test_utility [c] (n:i64) (newprices:[c]f64) (Ax:i64) : ?[ns][nd].trm.utility[ns][nd] =
  let [ns][nd] mp : trm.mp [n][c][Ax][ns][nd] = trm.mk n c Ax
  let mp = trm.set_newprices mp newprices
  let p = trm.simple_prices mp 0.85
  let tau = 1
  in trm.utility mp p tau

-- ==
-- entry: test_utility_2
-- input { 2i64 [100f64,100f64] 2i64  [[0.1f64, 0.1f64],[-0.1f64, -0.1f64]] }
-- output { [ [5.4f64, 4.500000002930908f64, 5.40000002528265f64, 4.500000002930908f64,
--             5.40000002528265f64, 8.50000015194252f64],
--            [f64.nan, -3.900000147521496f64, -3.000000125169754f64, -3.900000147521496f64,
--             -3.000000125169754f64, 0.10000000149011612f64],
--            [5.4f64, 4.500000002930908f64, 5.40000002528265f64, 4.500000002930908f64,
--             5.40000002528265f64, 8.50000015194252f64],
--            [f64.nan, -3.900000147521496f64, -3.000000125169754f64, -3.900000147521496f64,
--             -3.000000125169754f64, 0.10000000149011612f64],
--            [f64.nan, -4.000000149011612f64, -3.10000012665987f64, -4.000000149011612f64,
--             -3.10000012665987f64, 0.0f64]]
-- }

entry test_utility_2 [c] (n:i64) (newprices:[c]f64) (Ax:i64) (u_a_sq) : ?[ns][nd].trm.utility[ns][nd] =
  let [ns][nd] mp : trm.mp [n][c][Ax][ns][nd] = trm.mk n c Ax
  let mp = trm.set_newprices mp newprices
  let mp = trm.set_u_a_sq u_a_sq mp
  let p = trm.simple_prices mp 0.85
  let tau = 1
  in trm.utility mp p tau

-- ==
-- entry: test_bellman
-- input { 2i64 [100f64,100f64] 2i64 }
-- output { [8.670580092671319f64, 0.22769411217320684f64, 8.670580092671319f64,
--           0.22769411217320684f64, 0.12769411068309072f64] }

entry test_bellman [c] (n:i64) (newprices:[c]f64) (Ax:i64) : ?[ns].trm.ev[ns] =
  let [ns][nd] mp : trm.mp [n][c][Ax][ns][nd] = trm.mk n c Ax
  let mp = trm.set_newprices mp newprices
  let p = trm.simple_prices mp 0.85
  let tau = 1
  let util : trm.utility [ns][nd] = trm.utility mp p tau
  let tr = trm.age_transition mp
  let ev = trm.ev0 mp
  in trm.bellman mp util tr ev

-- ==
-- entry: test_bellmanN
-- input { 2i64 [100f64,100f64] 2i64 1i64 }
-- output { [8.670580092671319f64, 0.22769411217320684f64, 8.670580092671319f64,
--           0.22769411217320684f64, 0.12769411068309072f64] }
-- input { 2i64 [100f64,100f64] 2i64 2i64 }
-- output { [13.43965002284233f64, 5.039207394202223f64, 13.43965002284233f64,
--           5.039207394202223f64, 4.939207392712107f64] }
entry test_bellmanN [c] (n:i64) (newprices:[c]f64) (Ax:i64) (N:i64) : ?[ns].trm.ev[ns] =
  let [ns][nd] mp : trm.mp [n][c][Ax][ns][nd] = trm.mk n c Ax
  let mp = trm.set_newprices mp newprices
  let p = trm.simple_prices mp 0.85
  let tau = 1
  let util : trm.utility [ns][nd] = trm.utility mp p tau
  let tr = trm.age_transition mp
  in loop ev = trm.ev0 mp for _i < N do trm.bellman mp util tr ev
