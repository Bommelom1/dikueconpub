def bellman (x:[2]f64) : [2]f64 =
  let xs : *[2]f64 = replicate 2 0
  let xs[0] = x[1] ** 0.5
  let xs[1] = (1 - x[0] ** 2) ** 0.5
  in xs

def idd n i = tabulate n (f64.bool <-< (==i))

def run [n] (f: [n]f64 -> [n]f64) (x:[n]f64) : [n][n]f64 =
  tabulate n (jvp f x <-< idd n) |> transpose

def main (a:f64) : ([2]f64, [2][2]f64) =
  let x = replicate 2 a
  in (bellman x, run bellman x)
