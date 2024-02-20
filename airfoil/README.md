- NekRS Version, v23,1 (next)    
  https://github.com/Nek5000/nekRS/commit/7237c89d4ff232287f506b633963623efed2b12c

- From step 1, nek5000 and nekrs should get the same drag like this

  nekrs
  ```
    100 5.95619198e+00  1.09094709e-01  4.19992910e-02 6.70954183e-02 dragx(rs)
    100 5.95619198e+00  3.46084132e-01  3.42629908e-01 3.45422357e-03 dragy(rs)
    100 5.95619198e+00 -9.26109226e-02 -9.28473050e-02 2.36382423e-04 dragz(rs)
  ```

  nek5000 (compute in userhk)
  ```
     100  5.95619198461E+00  1.09094709316E-01  4.19992910125E-02  6.70954183038E-02      1dragx
     100  5.95619198461E+00  3.46084131851E-01  3.42629908282E-01  3.45422356810E-03      1dragy
     100  5.95619198461E+00 -9.26109225772E-02 -9.28473049997E-02  2.36382422531E-04      1dragz
  ```

- p55 is not supported
- variable viscosity is three but not tested

- TODO: merge to postprocessing
- one can uncomment the sij to print more. Nek5000's Sij has extra 2x.

- Case file is for naca airfoil case from Vishal. Mesh is not included here
