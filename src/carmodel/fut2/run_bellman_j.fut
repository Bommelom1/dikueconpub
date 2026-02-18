--
-- Solve the bellman equation for a given set of prices and parameters
--

import "autotrade/trmodel"
import "lib/github.com/diku-dk/linalg/dpsolve"

module mk_run (R:real) = {

  type t = R.t

  module trm = trmodel R
  module dps = mk_dpsolve_dense R

  def test (c:i64) (Ax:i64) : ([]t,[][]t) =
    let n = 1
    let tau = 0
    let newprices : [c]R.t = replicate c (R.i64 100)
    let [ns][nd] mp : trm.mp [n][c][Ax][ns][nd] = trm.mk n c Ax
    let mp = trm.set_newprices mp newprices
    let p = trm.simple_prices mp (R.f32 0.85)
    let util : trm.utility [ns][nd] = trm.utility mp p tau
    let tr = trm.age_transition mp
    let ev0 = trm.ev0 mp
    let f = trm.bellmanJ mp util tr
    in f ev0

  def run [c] (n:i64) (newprices:[c]t) (Ax:i64) (sa_max:i64) : ([]t,bool,i64,i64,i64,t) =
    let [ns][nd] mp : trm.mp [n][c][Ax][ns][nd] = trm.mk n c Ax
    let mp = trm.set_newprices mp newprices
    let p = trm.simple_prices mp (R.f32 0.85)
    let tau = 1
    let util : trm.utility [ns][nd] = trm.utility mp p tau
    let tr = trm.age_transition mp
    let ev0 = trm.ev0 mp
    let f = trm.bellmanJ mp util tr
    let param = dps.default with sa_max=sa_max
    let {res=ev,jac=_,conv=b,iter_sa=i,iter_nk=j,rtrips=k,tol} =
      dps.poly f ev0 param (R.f32 0)
    in (ev,b,i,j,k,tol)

  def runn [n][c][Ax][ns][nd] (mp:trm.mp[n][c][Ax][ns][nd]) (sa_max:i64) : ([n]t,bool,i64,i64,[n]i64) =
    let p = trm.simple_prices mp (R.f32 0.85)
    let tr = trm.age_transition mp
    let ev0 = trm.ev0 mp
    let param = dps.default with sa_max = sa_max
                            with sa_tol = R.f32 (1.0e-3)
    let comp (tau:i64) : (t,bool,i64,i64,i64) =
      let util : trm.utility [ns][nd] = trm.utility mp p tau
      let f = trm.bellmanJ mp util tr
      let {res=ev,jac=_,conv=b,iter_sa=i,iter_nk=j,rtrips=k,tol=_} =
	dps.poly f ev0 param (R.f32 0)
      let mx = R.maximum ev
      in (mx,b,i,j,k)
    let (mxs,bs,is,js,ks) =
      #[sequential_outer] map comp (iota n) |> unzip5
    in (mxs,
	reduce (&&) true bs,
	reduce (+) 0 is,
	reduce (+) 0 js,
	ks)
}

module R = f64

module r = mk_run R

entry main1 (c:i64) (Ax:i64) (sa_max:i64) : ([]R.t,bool,i64,i64,i64,R.t) =
  let pnews : [c]R.t = replicate c 100
  in r.run 2 pnews Ax sa_max

entry main (c:i64) (Ax:i64) (sa_max:i64) : []R.t =
  let pnews : [c]R.t = replicate c 100
  in (r.run 2 pnews Ax sa_max).0

entry mainn (n:i64) (c:i64) (Ax:i64) : ([n]R.t,bool,i64,i64,[n]i64) =
  #[unsafe]
  let sa_max = 3
  let pnews : [c]R.t = replicate c 100
  let u_0 = tabulate_2d n c
		(\i j -> R.(i64 5 + i64 2 * (i64 i + i64 j) / (i64 n + i64 c)))
  let mp = r.trm.mk n c Ax |> r.trm.set_u_0 u_0
  let mp = r.trm.set_newprices mp pnews
  in r.runn mp sa_max

entry test = r.test

entry tuneme (n:i64) : R.t =
  let c = 10
  let Ax = 15
  in mainn n c Ax |> (.0) |> reduce R.max (R.i64 0)

-- ==
-- entry: tuneme
-- input { 4i64 }
-- output { 166.996326f64 }

-- ==
-- entry: main
-- input { 2i64 2i64 200i64 }
-- output { [104.059715f64, 95.659248f64, 104.059715f64, 95.659248f64, 95.559250f64] }
-- input { 2i64 2i64 5i64 }
-- output { [104.059753f64, 95.659286f64, 104.059753f64, 95.659286f64, 95.559288f64] }

-- "echo 2 2 200 | ./run_bellman_j" converges, and so does
-- "echo 2 2 5 | ./run_bellman_j"!!
