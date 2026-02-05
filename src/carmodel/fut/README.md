
## Futhark Performance

We first compare running times for an implementation of the carmodel
using the successive approximation feature of the [`dpsolver`](util/dpsolver.fut) library. All measurements are performed on a
MacBook Pro (15-inch, 2016) with an 2,7 GHz Quad-Core Intel Core i7
processor, 16 GB of memory, and a Radeon Pro 460 4 GB graphics
card. Time values are measured as averages of 5 runs.

Comparing running times for a setup with four (4) households and
varying cartypes, each with a maxage of 20, we see that using Futhark
with the `opencl` and `multicore` backends run well compared to the
native `Matlab` version and the `c` Futhark backend:

![Comparing running times](/src/carmodel/fut/images/plot1.png)

(plot made using `make plot1.pdf`)

## Newton's Method

The [`dpsolver`](util/dpsolver.fut) library also allows for using
Newton's method for finding a fixpoint (with the aim of fast
convergence), as opposed to successively applying the Bellman
equation.

### Explicit Derivatives

Specifying a version of the Bellman function that returns not only the
successive element vector but also the derivative, in terms of a Jacobian
matrix, for the given element vector, leads to improved results:

![Comparing running times explicit Newton](/src/carmodel/fut/images/plot1j.png)

We are here using `dpsolver`'s *poly-solver*, which alternates between
using the successive approximation technique and Newton's method.  We
see that Matlab is quite efficient.

It turns out that we can tweak the Futhark implementation of the
LU-decomposition-based linear-system-solver by choosing a different
block-size than one. By choosing a block-size of 16, Futhark generates
improved code:

![Comparing running times explicit Newton Blocked 16](/src/carmodel/fut/images/plot1j_blocked16.png)

It would be great if we could make it up to auto-tuning to find the
optimal block-size (btw: some block-sizes - e.g., 32 - lead to wrong results.)

### Futhark AD

Unfortunately, using the `dpsolver`'s *poly_ad-solver*, which uses
Futhark's support for automatic differentiation in concert with using
the successive approximation technique and Newton's method leads to
poor performance:

![Comparing running times ad](/src/carmodel/fut/images/plot1ad.png)

The reason for this poor performance is slightly unclear. For the case
with four households, 20 car types, and a maximum car age of 25, the
Bellman function $F$, which is passed to the `dpsolver` fixpoint
solver, has type $R^{501} \rightarrow R^{501}$.  Using Newton's
method, LU-decomposition is perfomed a number of times on matrices of
size $501 \times 501$ (for solving a system of linear equations),
which in itself should not be a big problem (in fact, this part of the
computation works well with the non-ad approach). The problem seems to
be the use of automatic differentiation for computing the Jacobian
matrix (of size $501 \times 501$), which is done by applying Futhark's
`jvp` function to different direction vectors within a `tabulate`
construct:

```futhark
  def idd n i = tabulate n (\j -> if i==j then r.f32 1 else r.f32 0)

  def wrapj [n][m] (f: [n]t->[m]t) (x:[n]t) : ([m]t,[m][n]t) =
    (f x, #[sequential_outer] tabulate n (jvp f x <-< idd n) |> transpose)

  def poly_ad f x ap bet =
    let (y,_,b,i,j,k,tol) = poly (wrapj f) x ap bet
    in (y,b,i,j,k,tol)
```

One problem might be the lack of sufficient sharing (an excessive
amount of recomputation). For instance, the result of `f x` may be
computed many times, which may happen also for subexpressions). Notice
also the use of the `#[sequential_outer]` annotation, which was
necessary to avoid a drastic memory expansion.

The big question is how we (in general) can profile the code to find
out where time is spend and what memory is used for.

### Compiling and Running

To compile (c backend) and run with AD-generated Jacobian:

```
$ futhark c run_bellman_ad.fut
$ echo 4 15 25 | ./run_bellman_ad -e mainn -t /dev/stderr
3515810
[168.754394f64, 170.859001f64, 172.963673f64, 175.068405f64]
true
12i64
4i64
[1i64, 1i64, 1i64, 1i64]
```

To compile and run with hand-differentiated Jacobian:

```
$ futhark c run_bellman_j.fut
$ echo 4 15 25 | ./run_bellman_j -e mainn -t /dev/stderr
130156
[168.754394f64, 170.859001f64, 172.963673f64, 175.068405f64]
true
12i64
4i64
[1i64, 1i64, 1i64, 1i64]
```
