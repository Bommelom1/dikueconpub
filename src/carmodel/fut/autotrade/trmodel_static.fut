
import "trmodel"

module type trmodel_param = {
  val n  : i64
  val c  : i64
  val Ax : i64
}

module type trmodel_static = {

  type t

  val logsum [n] : (sigma:t) -> [n]t -> t

  -- model parameters
  include trmodel_param
  val ns : i64
  val nd : i64

  type mp
  val mp0 : mp

  type consumertype = i64

  -- car prices (for new and old cars)
  type prices = [Ax][c]t   -- c: ncartypes, Ax: maxage

  val simple_prices : mp -> prices

  -- utilities of transitions
  type utility = [ns][nd]t
  val utility : mp -> prices -> consumertype -> utility

  type~ transition -- = {trade:sp.mat[ns][ns], notrade:sp.mat[ns][ns]}
  val age_transition : mp -> transition

  -- expected values
  type ev = [ns]t
  val ev0 : mp -> ev

  -- the bellman equation
  val bellman : mp -> utility -> transition -> ev -> ev
}

module trmodel_static (R:real) (P:trmodel_param) : trmodel_static with t = R.t = {

  module B = trmodel R
  type t = B.t
  let logsum = B.logsum

  open P

  let mp0 [ns][nd] : B.mp [n][c][Ax][ns][nd] =
    let x = B.mk n c Ax
    in x

  type mp = B.mp [n][c][Ax][ns][nd]

  type consumertype = B.consumertype

  type prices = B.prices[c][Ax]

  let simple_prices = B.simple_prices

  type utility = B.utility[ns][nd]
  let utility = B.utility

  type~ transition = B.transition[ns]

  let age_transition = B.age_transition

  type ev = B.ev[ns]

  let ev0 = B.ev0

  let bellman = B.bellman
}
