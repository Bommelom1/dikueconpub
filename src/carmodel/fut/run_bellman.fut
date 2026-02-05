--
-- Solve the bellman equation for a given set of prices and parameters
--

import "autotrade/trmodel"
import "lib/github.com/diku-dk/linalg/dpsolve"

module mk_run (R:real) = {

  type t = R.t

  module trm = trmodel R
  module dps = mk_dpsolve_dense R

  def bellman = trm.bellman

  def run [c] (n:i64) (newprices:[c]t) (Ax:i64) (sa_max:i64) : ([]t,bool,i64,t) =
    let [ns][nd] mp : trm.mp [n][c][Ax][ns][nd] = trm.mk n c Ax
    let mp = trm.set_newprices mp newprices
    let p = trm.simple_prices mp (R.f32 0.85)
    let tau = 1
    let util : trm.utility [ns][nd] = trm.utility mp p tau
    let tr = trm.age_transition mp
    let ev0 = trm.ev0 mp
    let f = bellman mp util tr
    let param = dps.default with sa_max = sa_max
                            with sa_tol = R.f32 (1.0e-3)
    let {res=ev,conv=b,iter=i,tol,rtol=_} = dps.sa f ev0 param (R.f32 0)
    in (ev,b,i,tol)

  def runn [n][c][Ax][ns][nd] (mp:trm.mp[n][c][Ax][ns][nd]) (sa_max:i64) : ([n]t,bool,i64) =
    let p = trm.simple_prices mp (R.f32 0.85)
    let tr = trm.age_transition mp
    let ev0 = trm.ev0 mp
    let param = dps.default with sa_max = sa_max
                            with sa_tol = R.f32 (1.0e-3)
    let comp (tau:i64) : (t, bool, i64) =
      let util : trm.utility [ns][nd] = trm.utility mp p tau
      let f = bellman mp util tr
      let {res=ev,conv=b,iter=i,tol=_,rtol=_} = dps.sa f ev0 param (R.f32 0)
      let mx = R.maximum ev
      in (mx,b,i)
    let (mxs,bs,is) =
      #[sequential_outer] map comp (iota n) |> unzip3
    in (mxs,
	reduce (&&) true bs,
	reduce (+) 0 is)
}

module R = f64

module r = mk_run R

entry main1 (c:i64) (Ax:i64) (sa_max:i64) : ([]R.t,bool,i64,R.t) =
  let pnews : [c]R.t = replicate c 100
  in r.run 2 pnews Ax sa_max

entry main (c:i64) (Ax:i64) (sa_max:i64) : []R.t =
  let pnews : [c]R.t = replicate c 100
  in r.run 2 pnews Ax sa_max |> (.0)

entry mainn (n:i64) (c:i64) (Ax:i64) : ([n]R.t,bool,i64) =
  #[unsafe]
  let sa_max = 200
  let pnews : [c]R.t = replicate c 100
  let u_0 = tabulate_2d n c
		(\i j -> R.(i64 5 + i64 2 * (i64 i + i64 j) / (i64 n + i64 c)))
  let mp = r.trm.mk n c Ax |> r.trm.set_u_0 u_0
  let mp = r.trm.set_newprices mp pnews
  in r.runn mp sa_max

-- ==
-- entry: main
-- input { 2i64 2i64 200i64 }
-- output { [104.041557f64, 95.641090f64, 104.041557f64, 95.641090f64, 95.541092f64] }
-- input { 2i64 2i64 5i64 }
-- output { [26.364334f64, 17.963875f64, 26.364334f64, 17.963875f64, 17.863874f64] }

-- "echo 2 2 200 | ./run_bellman" converges, but not
-- "echo 2 2 5 | ./run_bellman"!!
