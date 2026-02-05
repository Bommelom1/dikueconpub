-- | Generic solver for Bellman equations. This module implements a
-- generic solver for Bellman equations in dynamic programming models
-- using John Rust's poly-algorithm (sa and newton are parts of this
-- module as well)
--
-- See Section 4 of https://editorialexpress.com/jrust/sdp/ndp.pdf
--
-- The module is enriched with a version of the solver that uses
-- Futhark's support for automatic differentiation (forward mode AD)
-- to compute the Jacobian for the passed Bellman function.

import "../lib/github.com/diku-dk/linalg/lu"

-- | Module type specifying generic solvers for Bellman equations.
module type dpsolver = {
  -- | The scalar type.
  type t

  -- | The type of solver parameters.
  type param = {
    sa_max          : i64,  -- Maximum and minimum numbers of
    sa_min          : i64,  --   successive approximation steps.
    sa_tol          : t,    -- Stopping tolerance for successive
                            --   approximation.
    max_fxpiter     : i64,  -- Maximum number of times to switch
                            --   between Newton-Kantorovich and
                            --   successive approximation iterations.
    pi_max          : i64,  -- Maximum number of Newton-Kantorovich steps.
    pi_tol          : t,    -- Final exit tolerance in fixed point
                            --   algorithm, measured in units of
                            --   numerical precision
    tol_ratio       : t     -- Relative tolerance before switching
                            --   to N-K algorithm when discount factor is
                            --   supplied as input in `poly`.
    }

  -- | The default parameter value
  val default : param

  -- | Type of structure
  type^  T 'a 'r
  val real : T t t
  val array : (n:i64) -> T ([n]t) ([n][n]t)
  val pair 'a 'b 'da 'db : T a da -> T b db -> T (a,b) ((da,db),(da,db))

  -- | Find a fixpoint for the function `f` using Newton-Kantorovich
  -- iterations with the initial guess `v`, and parameter `p`. The
  -- tolerance, the minimal, and the maximal number of iterations can
  -- be adjusted by altering the parameter `p`. The function `f`
  -- should return a pair of a new next approximation and the Jacobian
  -- matrix for the function `f` relative to the argument given. The
  -- function returns a quintuple containing an approximate fixpoint,
  -- a Jacobian matrix for the fixpoint, a boolean specifying whether
  -- the algorithm converged (according to the values in `p`), the
  -- number of iterations used, and finally, the tolerance of the last
  -- two fixpoint approximations (maximum of each dimension).
  val nk 'a 'd : T a d ->
		 (f: a -> (a,d)) -> (v:a) -> (p:param) ->
		 (a, d, bool, i64, t)

  -- | Find a fixpoint for the function `f` using a combination of
  -- successive approximation iterations and Newton-Kantorovich
  -- iterations. The initial guess is `v` and the parameter `p` is
  -- passed to the calls to `sa` and `nk`. The function `f` should
  -- return a pair of a new next approximation and the Jacobian matrix
  -- for the function `f` relative to the argument given. The function
  -- returns a 7-tuple containing an approximate fixpoint, a Jacobian
  -- matrix for the fixpoint, a boolean specifying whether the
  -- algorithm converged (according to the values in `p`), the number
  -- of iterations used for the total sa iterations, the total nk
  -- iterations, and the number of roundtrips. The 7'th element of the
  -- result tuple is the tolerance of the last two fixpoint
  -- approximations (maximum of each dimension).
  val poly 'a 'd : T a d -> (f: a -> (a,d)) -> (v:a) -> (p:param) ->
		 (b:t) -> (a, d, bool, i64, i64, i64, t)

  -- | Find a fixpoint for the function `f` using a combination of
  -- successive approximation iterations and Newton-Kantorovich
  -- iterations. The initial guess is `v` and the parameter `p` is
  -- passed to the calls to `sa` and `nk`. The function uses
  -- forward-mode automatic differentiation to compute the Jacobian
  -- matrix relative to the argument given to `f`. The function
  -- returns a 6-tuple containing an approximate fixpoint, a boolean
  -- specifying whether the algorithm converged (according to the
  -- values in `p`), the number of iterations used for the total sa
  -- iterations, the total nk iterations, and the number of
  -- roundtrips. The 6'th element of the result tuple is the tolerance
  -- of the last two fixpoint approximations (maximum of each
  -- dimension).
  val poly_ad 'a 'd : T a d -> (f:a->a) -> (v:a) -> (p:param) ->
		 (b:t) -> (a, bool, i64, i64, i64, t)
}

-- | Parameterised module for creating generic solvers.

module mk_dpsolver (r:real) : dpsolver with t = r.t = {

  type t = r.t
  local module lu = mk_lu r

  def blksz : i64 = 16  -- 1 or 16

  def ols [n] (m:[n][n]t) (b:[n]t) : [n]t =
    lu.ols blksz m b

  def eye (n:i64) (m:i64) : [n][m]t =
    tabulate_2d n m (\i j -> r.bool (i == j))

  def zero (n:i64) (m:i64) : [n][m]t =
    replicate n (replicate m (r.i64 0))

  type param = {
    sa_max          : i64,  -- Maximum number of contraction steps
    sa_min          : i64,  -- Minimum number of contraction steps
    sa_tol          : t,    -- Absolute tolerance before (in dpsolver.poly: tolerance before switching
                            --   to N-K algorithm)
    max_fxpiter     : i64,  -- Maximum number of times to switch between Newton-Kantorovich iterations
                            --   and contraction iterations.
    pi_max          : i64,  -- Maximum number of Newton-Kantorovich steps
    pi_tol          : t,    -- Final exit tolerance in fixed point algorithm, measured in units of
                            --   numerical precision
    tol_ratio       : t     -- Relative tolerance before switching to N-K algorithm
                            --   when discount factor is supplied as input in dpsolver.poly
  }

  def default : param =
    {sa_max      = 20,
     sa_min      = 2,
     sa_tol      = r.f64 1.0e-3,
     max_fxpiter = 35,
     pi_max      = 40,
     pi_tol      = r.f64 1.0e-13,
     tol_ratio   = r.f64 1.0e-03
    }


    module type typ = {
      type T 'a 'r
      val F64              : T f64 f64
      val Prod 'a 'r 'b 's : T a r -> T b s -> T (a,b) ((r,s),(s,r))
      val Arr        'a 'r : (n:i64) -> T a r -> T a ([n]r)

      val newton 'a 'd : T a d -> (a -> (a, d)) -> a -> a
    }

  module t : typ = {
    type T 'a 'r = {}
    def F64 = {}
    def Prod a b = {}
    def Arr n a = {}

    def newton t f x = x
  }

  def ex1 = t.(Prod (Arr 3 F64) (Arr 2 F64))

  -- | Type of structure
  type^  T 'a 'd =
    { eye : d
    , zero : d
    , mapd : (t -> t) -> d -> d
    , mapa : (t -> t) -> a -> a
    , maxa : a -> t
    , inv  : d -> d
    , mvm  : d -> a -> a
    , mmm  : d -> d -> d
    }

  def real : T t t =
    { eye = r.i64 1
    , zero = r.i64 0
    , mapd = \f x -> f x
    , mapa = \f x -> f x
    , maxa = \x -> x
    , inv = \x -> r.(i64 1 / x)
    , mvm = \a b -> r.(a * b)
    , mmm = \a b -> r.(a * b)
    }

  def array (n:i64) : T ([n]t) ([n][n]t) =
    { eye = tabulate_2d n n (\i j -> r.bool(i==j)
    , zero = replicate n (replicate n (r.i64 0))
    , mapd = \f -> map (map f)
    , mapa = map
    , maxa = r.maximum
    , }

  def pair 'a 'b 'da 'db : T a da -> T b db -> T (a,b) ((da,db),(da,db))


  def sa [m] (bellman : [m]t -> [m]t)
             (V0:[m]t)
             (ap:param)
             (bet:t) : ([m]t, bool, i64, t, t) =
      loop (V0,converged,i,tol,_rtol) = (V0, false, 0, r.i64 0, r.i64 0)
      while !converged && i < ap.sa_max do
        let V = bellman V0
	let tol' = reduce r.max (r.i64 0) (map2 (\a b -> r.(abs(a-b))) V V0)
	let rtol' = if i == 1 then r.i64 1
	 	    else r.(tol' / tol)
	let converged =
             -- Rule 1
             (i > ap.sa_min && (r.(abs(bet-rtol') < ap.tol_ratio)))
          || -- Rule 2
             --let adj = f64.(maximum V0 |> abs |> log10 |> ceil)
             --let ltol = ap.sa_tol * f64.(10 ** adj)
	     let ltol = ap.sa_tol
             in (i > ap.sa_min && r.(tol' < ltol))
	in (V, converged, i+1, tol',rtol')

  def nk [m] (bellman : [m]t -> ([m]t,[m][m]t))
             (V0:[m]t)
             (ap:param) : ([m]t, [m][m]t, bool, i64, t) =
    loop (V0,_dV0,converged,i,_tol) = (V0, zero m m, false, 0, r.i64 1)
      while !converged && i < ap.pi_max do
        let (V1, dV) = bellman V0
        --let V = map2 (-) V0 (la.matvecmul_row (la.inv dV) V1)
	let F = map2 (map2 (r.-)) (eye m m) dV
  	--let _hermitian =
	--  map2 (map2 (r.==)) F (transpose F) |> map (reduce (&&) true) |> reduce (&&) true
	let V = map2 (r.-) V0 (ols F (map2 (r.-) V0 V1))  -- NK-iteration
	-- do additional SA iteration for stability and accurate measure of error bound
	let (V0, _) = bellman V
	let tol' = r.maximum (map2 (\a b -> r.(abs(a-b))) V V0) -- tolerance

	-- adjusting the N-K tolerance to the magnitude of ev
	let adj = r.(maximum V0 |> abs |> log10 |> ceil)
	let ltol = r.(ap.pi_tol * (i64 10 ** adj)) -- Adjust final tolerance
        -- ltol=ap.pi_tol  -- tolerance

        let converged = r.(tol' < ltol) -- Convergence achieved
        in (V0, dV, converged, i+1, tol')

  -- dpsolver.poly: Solve for fixed point using a combination of Succesive Approximations (SA) and Newton-Kantorovich (NK) iterations
  --
  -- syntax:	(V, dV, iter) = dpsolver.poly(bellman, V0, ap, bet)
  --
  -- INPUT:
  --     bellman:   (V, P, dV) = Bellman equation with fixed point, V
  --	 	    Function on the form [V,dV]=bellman(V)
  --		    where V is the value function (m x 1), P is the policy function
  --		    and dV is the (m x m) Frechet derivative of the Bellman operator]
  --     V:  	    Initial guess value function,V.
  --		    [m x 1 matrix]
  --
  --     ap:	    Algorithm paramaters. See dpsolver.setup
  --
  --	 bet:	    Discout factor. Enters rule for stopping SA and switching to NK iterations.
  --		    SA should stop prematurely when relative tolerance is close to bet.
  --
  -- OUTPUT:
  --     V:         m x 1 matrix. Fixed point, V
  -- 	 P:	    Policy function at fixed point
  -- 	 dV:	    Frechet derivative of the Bellman operator

  def poly [m] (bellman : [m]t -> ([m]t,[m][m]t))
               (V0:[m]t)
               (ap:param)
               (bet:t) : ([m]t,[m][m]t,bool,i64,i64,i64,t) =
    loop (V0,_dV,converged,i,j,k,_tol) = (V0, zero m m, false, 0, 0, 0, r.i64 1)
      while !converged && k < ap.max_fxpiter do
        -- poly-algorithm loop (switching between sa and nk)
        let (V1,_,i',_,_) = sa ((.0) <-< bellman) V0 ap bet
	let (V2, dV, c2, j', tol) = nk bellman V1 ap
        in (V2,dV,c2,i+i',j+j',k+1,tol)


  def idd n i = tabulate n (\j -> if i==j then r.f32 1 else r.f32 0)

  def wrapj [n][m] (f: [n]t->[m]t) (x:[n]t) : ([m]t,[m][n]t) =
    (f x, #[sequential_outer] tabulate n (jvp f x <-< idd n) |> transpose)

  def poly_ad f x ap bet =
    let (y,_,b,i,j,k,tol) = poly (wrapj f) x ap bet
    in (y,b,i,j,k,tol)

}
