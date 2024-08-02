# Computation of Nusselt number
This compute the averaged Nusselt number on given surface(s)

### Formulation

$$Nu = \frac{h L}{k}$$

- $L$: Characteristic length
- $k$: Thermal conductivity
- $h$: Heat transfer coefficient

$$h = \frac{-k_{fluid} \nabla\T\cdot{n}}{T_s-T_\infty}$$

Under non-dimensionalization, $\hat{T} = \frac{T-T_{\infty}}{\delta T}$ with $\deltaT={T_s-T_\infty}$, the averaged Nusselt number is computed as the surface integral

$$\bar{Nu} = \frac{1}{|S|} \int_S \nabla T \cdot \vec{n} dS$$
where $S$ is the target surface and $|S|$ is its area.

### Verson

- NekRS version: v23

### Usage
This is the turbPipe example with unit heat flux on the surface.

- Include the file in udf
  ```
  #include "comp_nusselt.hpp"
  ```

- Build kernel in `UDF_LoadKernels`
  ```
  void UDF_LoadKernels(occa::properties &kernelInfo)
  {
    nusselt::buildKernel(kernelInfo);
  }
  ```

- Allocate array size in `UDF_Setup`
  ```
  nusselt::setup(nrs->meshV, nrs->fieldOffset);
  ```

- Compute and print with the number with the following call
  ```
  bool print = true;                // print the number 
  std::vector<int> bidWall = {3};   // boundary id for surface(s)
  auto nu = nusselt::compute(nrs, bidWall, nrs->cds->o_S, time, tstep, print);
  ```


### Verification


It matches with the Nek5000 number (see `turbPipe.usr/userchk`)
```
     step        time    integral        area          Nu
      200  1.2000e+00  3.1407e+05  6.2832e+01  4.9985e+03 nusselt(rs)
      200  1.2000E+00  3.1407E+05  6.2832E+01  4.9985E+03 nusselt(5k)
```

TODO: compare is with reference (DNS, expriments, and correlation)?

