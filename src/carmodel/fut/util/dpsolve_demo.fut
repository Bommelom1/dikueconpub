-- ---
-- title: Dynamic Programming in Futhark
-- author: Martin Elsman
-- date: October 23, 2021
-- citation-style: acm-sig-proceedings-long-author-list
-- ...

-- This note describes how the Futhark [@futhark2017;
-- @futharkbook2018] dynamic programming library `dpsolve`[^1] can be
-- used to find fix-points for functions from $R^n$ to $R^n$ and how
-- solutions to multiple instances of a dynamic programming problem
-- can be computed in parallel on GPUs. We give both
-- single-dimensional and multi-dimensional examples and we show how
-- the Futhark automatic differentiation feature may relieve
-- programmers from specifying explicitly the Jacobian matrices, which
-- are necessary for using `dpsolve`'s fast converging Newtonian
-- functionality.

-- [^1]: The Futhark library `dpsolve` is based on a Matlab library
-- implemented by Bertel Schjerning, ECON, University of Copenhagen.

-- # Introduction
--
-- A standard approach for finding fix-points for numerical functions
-- from $R^n$ to $R^n$ is to use the technique of successive
-- approximations. Following Section 4 of [Numerical Dynamic
-- Programming in
-- Economics](https://editorialexpress.com/jrust/sdp/ndp.pdf), by John
-- Rust [@RUST1996619], the dynamic programming solver that we shall
-- apply here first uses a number of successive approximation steps
-- before it applies a more efficient Newtonian method for narrowing
-- in on a fix-point. The latter method requires that the user
-- specifies how to compute the Jacobian matrix (of type $R^{n \times
-- n}$) given an approximate fix-point. The Jacobian is then computed
-- for each Newtonian step.

-- # Example: Intersection of a circle and a quadratic equation
--
-- Following the example in [Jim Lambers' MAT 461/561 lecture
-- notes](https://www.scribd.com/document/685401140/Lecture-22), we first set
-- out to find the intersection between the unit circle ($x_1^2+x_2^2 = 1$) and
-- the quadratic equation $x_2 = x_1^2$.[^2]
--
-- [^2]: Analytically, the solution can easily be found by solving the
-- quadratic equation $x_2^2 + x_2 - 1 = 0$, which leads to the
-- solution $x_1 = 0.786151377757$ and $x_2 = 0.61803398874989$.
--
-- We first define an operator for which we want to find a
-- fix-point. To ensure that the natural matrix norm of the Jacobian
-- matrix for the function is less than 1 ($0 \leq x_1 \leq 1$ and $0
-- \leq x_2 \leq 1$), we give the following definition of the fix-point
-- operator $G$:
-- \begin{align} G (x_1,x_2) &= (\sqrt{x_2}, \sqrt{1-x_1^2}) \label{eq:1}
-- \end{align}
-- Without having to define the Jacobian matrix for the function, we
-- can find an approximation to the fix-point using the successive
-- approximation functionality of the `dpsolve` library.

-- We first `import` the library `dpsolve` and instantiate the parameterised
-- module `mk_dpsolve_dense` to the `f64` representation of floats:

import "../lib/github.com/diku-dk/linalg/dpsolve"
module dps = mk_dpsolve_dense f64

-- The function `dps.sa` that we shall apply has the following type:
-- ```futhark
-- val sa [m] : (f:[m]t->[m]t) -> (v:[m]t) -> (p:param) -> (b:t)
--              -> {res:[m]t,conv:bool,iter:i64,tol:t,rtol:t}
-- ```
--
-- Here `t` is identical to `f64` due to the `f64` module
-- instantiation of the `mk_dpsolve_dense` parameterised module.

-- We now define the `bellman` equation for which we want to find a
-- fix-point:

def bellman (x:[2]f64) : [2]f64 =
  [f64.sqrt x[1], f64.sqrt(1 - x[0] ** 2)]

-- Notice that we use projections from the argument vector to access the scalar
-- values.  The following Futhark function makes a call to the `dps.sa` function
-- with the above `bellman` function given as a parameter:

def test_sa (sa_max:i64) (sa_tol:f64) : ([2]f64, bool, i64, f64) =
  let v0 = [0.5, 0.5]
  let ap = dps.default with sa_max = sa_max
                       with sa_tol = sa_tol
  let {res, conv=b, iter=i, tol, rtol=_} = dps.sa bellman v0 ap 0
  in (res, b, i, tol)

-- The function `dps.sa` also takes an initial guess as argument (`v0`) together
-- with an `ap` value that defines some slightly modified default parameter
-- settings (max iterations, max tolerance, etc.)
--
-- We can now call the function:

-- > test_sa 60i64 1e-3

-- We see that after 56 iterations, a fix-point is found with a
-- tolerance below `1e-3`, meaning that the last iteration step
-- contributed to a change in value of less than `1e-3` for both $x_1$
-- and $x_2$. For improved precision, many more iterations are
-- required:

-- > test_sa 200i64 1e-9

-- # Faster convergence with Newton's method
--
-- The function $G$, as defined in $\eqref{eq:1}$, has the following Jacobian matrix:
-- \begin{align*} J_G(x_1,x_2) &= \left [ \begin{array}{cc}
--                              0                           & \frac{1}{2\sqrt{x_2}} \\
--                              \frac{-x_1}{\sqrt{1-x_1^2}} & 0
--                            \end{array} \right
--                         ]
-- \end{align*}

-- The following version of the `bellman` function takes its input as
-- an array of size 2 and returns, along with the function result, the
-- Jacobian matrix, relative to the argument:

def bellman_j (a:[2]f64) : ([2]f64, [2][2]f64) =
  let x1 = a[0]
  let x2 = a[1]
  let res = [f64.sqrt x2, f64.sqrt(1-x1**2)]
  let j = [[0                       , 1/(2*f64.sqrt x2) ],
           [-x1/(f64.sqrt(1-x1**2)) , 0                 ]]
  in (res, j)

-- The function `dps.poly` that we shall apply has the following type:
-- ```futhark
-- val poly [m]: (f: [m]t -> ([m]t,[m][m]t)) -> (v:[m]t) -> (p:param) -> (b:t)
--               -> {res:[m]t,jac:[m][m]t,conv:bool,iter_sa:i64,iter_nk:i64,rtrips:i64,tol}
-- ```
--
-- Again, here `t` is identical to `f64` due to the `f64` module instantiation
-- of the `mk_dpsolve` parameterised module. The function finds a fix-point for
-- the function `f` using a combination of successive approximation iterations
-- and Newton-Kantorovich iterations. The initial guess is `v` and the parameter
-- `p` controls the iteration passes. The function `f` should return a pair of a
-- new next approximation and the Jacobian matrix for the function `f` relative
-- to the argument given. The function returns a record containing an
-- approximate fix-point, a Jacobian matrix for the fix-point, a boolean
-- specifying whether the algorithm converged (according to the values in `p`),
-- the total number of iterations used for `sa` iterations, the total number of
-- Newton-Kantoovich iterations, the number of round-trips, and the tolerance of
-- the last two fix-point approximations (maximum of each dimension).

def test_poly (sa_max:i64) : ([2]f64, bool, i64, i64, i64, f64) =
  let v0 = [0.5, 0.5]
  let ap = dps.default with sa_max = sa_max
  let {res, jac=_, conv, iter_sa, iter_nk, rtrips, tol} = dps.poly bellman_j v0 ap 0
  in (res, conv, iter_sa, iter_nk, rtrips, tol)

-- > test_poly 5i64

-- Notice that the programmer has manually provided code for computing the
-- Jacobian matrix for the function. The result is a fix-point with a tolerance
-- below 1e-15, computed with an initial number of 5 successive approximation
-- iterations followed by 4 Newtonian iterations (1 round-trip was used).

-- # Futhark AD

-- We can relieve the programmer from manually providing the code for the
-- Jacobian matrix by using the automatic differentiation feature of Futhark,
-- which provides a function
-- [`jvp`](https://github.com/diku-dk/futhark/issues/1249) that performs
-- forward-mode automatic differentiation on arbitrary Futhark functions. An
-- alternative is to encode float computations using so-called dual-numbers,
-- following the approach of [AD with dual
-- numbers](https://futhark-lang.org/examples/dual-numbers.html), but we shall
-- not dive into this possibility here.
--
-- We first define a function `wrapj` that takes a function of type
-- `[n][m].[n]f64->[m]f64` and turns it into a function of type
-- `[n][m].[n]f64->([m]f64,[m][n]f64)` that, besides from the function result,
-- returns the Jacobian matrix of the function:

def idd n i = tabulate n (f64.bool <-< (==i))

def wrapj [n][m] (f: [n]f64->[m]f64) (x:[n]f64) : ([m]f64,[m][n]f64) =
  (f x, tabulate n (jvp f x <-< idd n) |> transpose)

-- Functions wrapped with the `wrapj` function can now be used directly with the
-- `dps.poly` function. Let's try it out in practice:

def test_poly_jvp (sa_max:i64) : ([2]f64,bool,i64,i64,i64,f64) =
  let v0 = [0.5, 0.5]
  let ap = dps.default with sa_max = sa_max
  let {res, jac=_, conv, iter_sa, iter_nk, rtrips, tol} =
    dps.poly (wrapj bellman) v0 ap 0
  in (res, conv, iter_sa, iter_nk, rtrips, tol)

-- > test_poly_jvp 5i64

-- We see that we get the same results with `test_poly_jvp` as we get with
-- `test_poly`.

-- To make it even easier for the programmer, the `dps` module includes a
-- version of the `poly` function, called `polyad`, that takes care of computing
-- the Jacobian of the passed function:

def test_polyad (sa_max:i64) : ([2]f64,bool,i64,i64,i64,f64) =
  let v0 = [0.5, 0.5]
  let ap = dps.default with sa_max = sa_max
  let {res, conv, iter_sa, iter_nk, rtrips, tol} = dps.polyad bellman v0 ap 0
  in (res, conv, iter_sa, iter_nk, rtrips, tol)

-- > test_polyad 5i64

-- # Going parallel
--
-- The iterative approaches that the `dpsolve` functionality implements for
-- finding fix-points are inherently sequential, except from the matrix
-- operations applied in the Newton-Kantorovich iterations (assuming a
-- high-number of dimensions). Instead of parallelising the actual fix-point
-- resolution, we shall see how we can find many fix-points in parallel, which
-- is sometimes a useful approach for speeding up an application.
--
-- Following up on the task of finding intersection points between a circle and
-- a simple quadratic equation, let us investigate how the x-dimension of the
-- intersection points changes when the circle radius increases.
--
-- We first parameterise the bellman equation over the radius of the circle:

def bellmanr (r:f64) (a:[2]f64) : ([2]f64, [2][2]f64) =
  let f (a:[2]f64) = [f64.sqrt a[1], f64.sqrt(r**2-a[0]**2)]
  let res = f a
  let j = [[0                            , 1/(2*f64.sqrt a[1]) ],
           [-a[0]/f64.sqrt(r**2-a[0]**2) , 0                   ]]
  in (res, j)

-- We then create a function that implements an outer map over a call to
-- `dps.poly` with varying radius:

def linspace (n: i64) (start: f64) (end: f64) : [n]f64 =
  tabulate n (\i -> start + f64.i64 i * ((end-start)/f64.i64 n))

def test_polyr (n:i64) (sa_max:i64) : (bool, i64, [n]f64, [n]f64) =
  let ap = dps.default with sa_max = sa_max
  let rs = linspace n 1 20
  let ress = map (\r -> let v0 = [0.5,0.5]
                        let {res, jac=_, conv, iter_sa, iter_nk, rtrips=_, tol=_} =
                          dps.poly (bellmanr r) v0 ap 0
                        in (res,res[0],conv,iter_sa+iter_nk)) rs
  let converged = reduce (&&) true (map (.2) ress)
  let xs = map (.1) ress
  let iterations = reduce (+) 0 (map (.3) ress)
  in (converged, iterations, rs, xs)

-- Here is a call to `test_polyr` with 4 different radius values (between 1 and
-- 20) and an `sa_max` value of 3:

-- > test_polyr 4i64 3i64

-- We can use the plot functionality of [literate
-- Futhark](https://futhark-lang.org/examples/literate-basics.html) to plot 1000
-- points relating radius values with associated found $x$-values (and compare
-- it with a plot of the `sqrt`-function):

def test_polyr_rxs (n:i64) (sa_max:i64) : ([n]f64, [n]f64) =
  test_polyr n sa_max |> (\(_,_,rs,xs) -> (rs,xs))

def xys f n start end =
  unzip (map (\x -> (x, f x)) (linspace n start end))

def sqrt_coords = xys f64.sqrt

-- > :plot2d {rxs=test_polyr_rxs 1000i64 3i64,
--            sqrt=sqrt_coords 1000i64 1f64 21f64};
-- size:(1000, 400)

-- # A few single-dimensional examples

-- We now consider a single-dimensional case, for which we want to find the $x$
-- for which $f(x) = \cos x$.

def test_poly1d (sa_max : i64) : ([1]f64, [1][1]f64, bool, i64, i64, i64, f64) =
  let ap = dps.default with sa_max = sa_max
  let {res, jac, conv, iter_sa, iter_nk, rtrips, tol} =
    dps.poly (\x -> ([f64.cos x[0]],
                     [[- f64.sin x[0]]]))
             [0.7] ap 0
  in (res, jac, conv, iter_sa, iter_nk, rtrips, tol)

-- > test_poly1d 0i64

-- For another example, we want to compute $\sqrt{2}$ by finding the
-- fix-point to the equation $f(x)=\frac{1}{2}(x+\frac{2}{x})$.

def test_sqrt (sa_max : i64) : ([1]f64, [1][1]f64, bool, i64, i64, i64, f64) =
  let ap = dps.default with sa_max = sa_max
  let {res, jac, conv, iter_sa, iter_nk, rtrips, tol} =
    dps.poly (\x -> ( [ 0.5 * (x[0]+2/x[0]) ],
                      [[ 2*x[0] ]] )
             ) [1.4] ap 0
  in (res, jac, conv, iter_sa, iter_nk, rtrips, tol)

-- > test_sqrt 0i64

-- Remarkably, in 4 steps we reach a fix-point of 1.41421356237... with
-- a tolerance of 1.11e-15.

-- # Solving Systems of Linear Equations

-- It is common to use iterative methods also to solve systems of linear
-- equations. Such methods include the Jacobi method and the Gauss-Seidel method
-- [@IterativeMethods2003, p.127]. In general, a system of linear equations can
-- be written on the form $A\mathbf{x} = \mathbf{b}$, where $A$ is a known
-- square matrix, $\mathbf{b}$ is a known vector, and $\mathbf{x}$ is the vector
-- we seek to find. Provided we split up $A$ into a lower-triangular matrix $L$,
-- an upper-triangular matrix $U$, and a diagonal matrix $D$, such that
-- $A=L+U+D$, it turns out that we can write a recurrence equation for the
-- linear system on the form \begin{align} \mathbf{x}_{k+1} &= G\mathbf{x}_k +
-- \mathbf{f} \label{eq:2} \end{align} where, for the Jacobi method, we further
-- have $G = I - D^{-1}A$ and $\mathbf{f}=D^{-1}\mathbf{b}$. Here is Futhark
-- code that implements the Jacobi method:

def dotprod [n] (u:[n]f64) (v:[n]f64) : f64 =
  reduce (+) 0.0 (map2 (*) u v)

def matvecmul [n][m] (A: [n][m]f64) (v: [m]f64) : [n]f64 =
  map (dotprod v) A

def matmul [n][p][m] (us: [n][p]f64) (vs: [p][m]f64) : [n][m]f64 =
  map (matvecmul (transpose vs)) us

def binop [n] (f:f64->f64->f64) (a:[n][n]f64) (b:[n][n]f64) : [n][n]f64 =
  map2 (map2 f) a b

def diag_ex [n] (A:[n][n]f64) : [n]f64 =
  map (\i -> A[i][i]) (iota n)

def diag [n] (a:[n]f64) : [n][n]f64 =
  tabulate_2d n n (\i j -> if i == j then a[i] else 0.0)

def jacobi [n] (A:[n][n]f64) (b:[n]f64) : [n]f64 -> [n]f64 =
  let D' = diag_ex A |> map (1.0/) |> diag -- D^{-1}
  let I = diag (tabulate n (\_ -> 1.0))
  let G = binop (-) I (matmul D' A)
  let f = map (\i -> D'[i][i]*b[i]) (iota n)
  in \x -> map2 (+) (matvecmul G x) f

-- We can try out the Jacoby method by attempting to find a solution to the
-- following system of linear equations:
-- \begin{align} 3x_0 +  x_1  -  x_2  &= 1 \label{eq:3} \\
--                x_0 -  x_1  + 2x_2  &= 8 \nonumber \\
--               -x_0 +  x_1  - 3x_2  &= -1 \nonumber
-- \end{align}
-- Here is a function that iterates the Jacobi method a number of times
-- starting with an initial constant vector:
def test_jacobi [n] (A:[n][n]f64) (b:[n]f64) (k:i64) : [n]f64 =
  let pow f k x = loop x for _i < k do f x
  let f = jacobi A b
  in pow f k (replicate n 1f64)

def A = [[3f64,1,-1],[1.0,-1,2],[-1.0,1,-3]]
def b = [1f64,8,-1]

-- > test_jacobi A b 5

-- With a few more iterations, we converge towards a desired solution:

-- > test_jacobi A b 50

-- > test_jacobi A b 51

-- We can test that the found solution is close to be correct:

def test_jacobi_ok : [3]f64 = matvecmul A (test_jacobi A b 51)

-- > test_jacobi_ok

-- Now, let's try instead to use the function `polyad` to find a solution:

def test_jacobi_polyad (sa_max: i64) : ([3]f64,bool,i64,i64,i64,f64) =
  let f = jacobi (copy A) (copy b)
  let x0 = replicate 3 1f64
  let ap = dps.default with sa_max = sa_max
  let {res, conv, iter_sa, iter_nk, rtrips, tol} = dps.polyad f x0 ap 0
  in (res, conv, iter_sa, iter_nk, rtrips, tol)

-- > test_jacobi_polyad 2i64

-- In two steps we reach a fix-point with a tolerance of 3.5e-15. We shall not
-- here go into details about under what conditions the recurrence converges to
-- a fix-point. For a proper analysis of these aspects, we refer the reader to
-- [@IterativeMethods2003], which also presents improvements to the Jacobi
-- method in terms of the Gauss-Seidel variation and the successive over
-- relaxation (SOR) method. An interesting aspect of these methods is that they
-- may work well also for very large but sparse matrices. In particular, we
-- notice that, for the Jacobi method, the matrix $G$ has similar sparsity
-- structure as $A$, provided $A$ is non-singular (i.e., have non-zero diagonal
-- elements). However, as fun as this exercise may be, notice that, for each
-- Newton-Kantorovich iteration, the `polyad` function will solve a linear
-- equation system (using the `ols` function from the `linalg` module), which is
-- as big as the equation system that we aim at solving, thus, not much is
-- gained in the case of solving linear equation systems. It is really for the
-- case of solving non-linear systems of equations that the technique is
-- valuable.

-- # Conclusion

-- We have seen how we can use the `dpsolve` library to solve
-- multi-dimensional fix-point equations. We have also seen how we can
-- solve multiple problems in parallel using Futhark's second-order
-- array combinators.

-- # References
