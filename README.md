# dikueconpub

Public code base for solving economic problems

- `src/logit`: MC code in Python and Futhark for computing some economic measures. To run the Futhark code:
  ```
  $ cd src/logit
  $ make
  ...
  ```

- `src/carmodel`: Matlab and Futhark code for finding the equilibrium of a
  used-car model, using dynamic programming. For a description of the problem
  see [Equilibrium Trade in Automobiles](https://github.com/melsman/dikueconpub/blob/main/src/carmodel/docs/eqb.pdf). To run the Futhark code:
  ```
  $ cd src/carmodel/fut
  $ make
  ...
  $ make plots
  ...
  ```

  On an Apple M2 Max Macbook Pro (2023), the generated [plots](https://github.com/melsman/dikueconpub/blob/main/src/carmodel/fut/plots-mac.pdf) show promising results for Futhark...
