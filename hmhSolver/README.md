# Example of solving extra Helmholtz solver

- version: v24 (repo/next 07/10/24 + some fixes that will be released later)

$$-\nabla\cdot(h_1\nabla T) + h_2 T = rhs$$

### Settings

We use Scalar01, Scalar02, Scalar03 to solve extra Helmholtz solvers.
The Scalar00 (Temperature) solves the same things as the Robin example.

```
s01: Dirichlet BC
s02: Neumann BC
s03: Robin BC
```

Each tests are solved twice where at
```
step= 1: solve Poisson equation (h2=0)
step= 11: solve helmholtz equation
```


### API

To use the slot of `scalarId`, put a session in par and set tol, filter, solver, and preconditoner there.

```
hmhSolver::setup(nrs, scalarId);
   Some checks + initialize the timer. Most part is done by NekRS default. There are not much to setup actually. 
   After calling setup, one can set cds->compute[scalarId] = 0 to turn off scalar solver

hmhSolver::solve(nrs, scalarId, time, tstep, o_h1, o_h2, o_rhs);
   This is the main driver, one can enter arbitray (h1,h2,rhs)
   The BC will use the scalar's oudf + implicit linear term
```


### Verification (WIP)

```
$ grep ^hmh_test logfile |tail
hmh_test poi  step=1 time=0.005 scalarId=1 (h1,h2)=(0.101321,0) PGMRES iter=63 err(abs,rel)=(4.48281e-08,2.41527e-09) res=(3.00e+02,6.32e-15)
hmh_test poi  step=1 time=0.005 scalarId=2 (h1,h2)=(0.101321,0) PGMRES iter=119 err(abs,rel)=(8.92807e-08,2.41493e-09) res=(1.84e+02,9.65e-15)
hmh_test poi  step=1 time=0.005 scalarId=3 (h1,h2)=(0.101321,0) PGMRES iter=90 err(abs,rel)=(4.48798e-08,1.12007e-09) res=(2.41e+02,7.40e-15)
hmh_test hmh  step=11 time=0.055 scalarId=1 (h1,h2)=(0.101321,-4) PGMRES iter=39 err(abs,rel)=(4.48253e-08,4.83319e-09) res=(8.50e-09,8.55e-15)
hmh_test hmh  step=11 time=0.055 scalarId=2 (h1,h2)=(0.101321,-4) PGMRES iter=500 err(abs,rel)=(8.96332e-08,2.42146e-09) res=(8.41e-09,1.36e-13)
hmh_test hmh  step=11 time=0.055 scalarId=3 (h1,h2)=(0.101321,-4) PGMRES iter=95 err(abs,rel)=(4.4811e-08,1.35542e-09) res=(3.41e-09,9.82e-15)
```
