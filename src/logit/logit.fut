--- First some utility/boilerplate.

import "lib/github.com/diku-dk/cpprandom/random"
module rng = xorshift128plus
module rand_f64 = normal_distribution f64 rng
type rng = rng.rng

let rand_vector (r: rng) (N: i64) =
  let (rngs, vs) = r
                   |> rng.split_rng N
                   |> map (rand_f64.rand {mean=0, stddev=1})
                   |> unzip
  in (rng.join_rng rngs, vs)

let rand_matrix (r: rng) (N: i64) (K: i64) =
  let (rngs, vs) = r
                   |> rng.split_rng (N*K)
                   |> map (rand_f64.rand {mean=0, stddev=1})
                   |> unzip
  in (rng.join_rng rngs, unflatten vs)

let dotprod xs ys = map2 (*) xs ys |> f64.sum

let matmul [n][m][k] (xss: [n][m]f64) (yss: [m][k]f64) =
  map (\xs -> map (dotprod xs) (transpose yss)) xss

let mean [n] (xs: [n]f64): f64 =
  map (/f64.i64 n) xs |> f64.sum

--- Now a more or less straightforward translation of the NumPy code.

let compute_utils [N][Kx][Kz][J]
                  (x: [N][Kx]f64)
                  (z: [Kz][J]f64)
                  (beta: [Kx][Kz]f64)
                  (MAX_RESCALE: bool)
                  : [N][J]f64 =
  let U = x `matmul` beta `matmul` z
  in if MAX_RESCALE then let Umax = map f64.maximum U
                         in map2 (\Umax' -> map (\x -> x - Umax')) Umax U
     else U

let ccps [N][Kx][Kz][J]
         (x: [N][Kx]f64)
         (z: [Kz][J]f64)
         (beta: [Kx][Kz]f64) =
  let EU = map (map f64.exp) (compute_utils x z beta true)
  let prob = map (\row -> map (/f64.sum row) row) EU
  in prob

let simulate_choices [N][J] (rng: rng)
                            (ccp: [N][J]f64)
                          : [N]i64 =
  let ccp_cum = map (scan (+) 0) ccp
  let (_, u) = rand_vector rng N
  in loop d = replicate N 0 for j < J-1 do
     let I = map3 (\i x y -> if x > y then i else -1)
                  (iota N) u ccp_cum[:,j]
     in scatter d I (replicate N (j+1))

let loglikelihood [N][Kx][Kz][J]
                  (x: [N][Kx]f64)
                  (z: [Kz][J]f64)
                  (d: [N]i64)
                  (beta: [Kx][Kz]f64)
                 : [N]f64 =
  let U = compute_utils x z beta true
  let logsum = map (map f64.exp >-> f64.sum >-> f64.log) U
  let U_d = map2 (\Ur i -> Ur[i]) U d
  in map2 (-) U_d logsum

let neg_ll [N][Kx][Kz][J]
           (beta: []f64)
           (x: [N][Kx]f64)
           (z: [Kz][J]f64)
           (d: [N]i64)
         : f64 =
  let B = unflatten (beta:[Kx*Kz]f64)
  in -mean(loglikelihood x z d B)

let simulate_demographics (rng: rng) (N: i64) (J: i64) (Kz: i64) (Kx: i64) =
  let x = replicate Kx 1 |> replicate N
  let (rng, xs) = rand_matrix rng N (Kx-1)
  let x[:,1:] = xs
  let (rng, z) = rand_matrix rng Kz J
  in (rng, x, z)

--- Data generation/entry points.

-- Use to generate a data file, e.g:
--
--   $ futhark opencl logit.fut
--   $ echo 123 10000 3 3 2 | ./logit -e generate_data -b > logit.in
--   $ echo 123 10000 3 3 2 | ./logit -d AMD -e generate_data -b > logit.in
entry generate_data (seed: i32) (N: i64) (J: i64) (Kz: i64) (Kx: i64) =
  let rng = rng.rng_from_seed [seed]
  let (rng, x, z) = simulate_demographics rng N J Kz Kx
  let (rng, beta) = rand_matrix rng Kx Kz
  let ccp = ccps x z beta
  let d = simulate_choices rng ccp
  in (x, z, beta, ccp, d)

-- Use to time the criterion function, e.g:
--
--  $ ./logit -e criterion < logit.in  -t /dev/stderr -r 10
--  $ ./logit -d AMD -e criterion < logit.in  -t /dev/stderr -r 10
entry criterion (x: [][]f64)
                (z: [][]f64)
                (beta: [][]f64)
                (_ccp: [][]f64)
                (d: []i64) =
  neg_ll (flatten beta) x z d
