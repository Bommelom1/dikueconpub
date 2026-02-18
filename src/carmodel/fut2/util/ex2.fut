def bellman2 (x:[2]f64) : ([2]f64, [2][2]f64) =
  let a = x[1] ** 0.5
  let b = (1 - x[0] ** 2) ** 0.5
  in ([a, b],
      [[0, 1/(2*a)],
       [-x[0]/b, 0]]
      )

def main (a:f64) : ([2]f64, [2][2]f64) =
  bellman2 [a,a]
