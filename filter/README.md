# Explicit filter

Here we apply the "explicit" filter used in Nek5000.
The case is built on the top of the either example.
TODO: maybe change to an example that does need filtering.

- Version: NekRS v24-pre (repo/next, tested on 08/18/24)


### Formulation

See the Nek5000 [docuentation](https://nek5000.github.io/NekDoc/problem_setup/filter.html).

Basically, in each direction, it changes the basis, rescale the magnitude with diagonal matrix then change the basis back. The whole action is stored intoa 1D matrix and the kernel simply do max-vec in each direction, elements by elements.


### Case setup

1. Put the files `nekFilter.hpp` and `nekFilter.okl` to the case folder.

2. `#include "nekFilter.hpp"` in UDF file.

3. Turn off `regularization` for the desired field.       
   This is user's responsibility.
   We don't check this in the code since user can have both filter turning on.

4. Call `nekFilter::buildKernel(kernelInfo);` in the UDF function `UDF_LoadKernels`.

5. Call the setup function with arbitary "tag", following by the mesh, filterNode and filterWeight. It can have multiple filter for each fields.
   ```
   nekFilter::setup("vel", nrs->meshV, 1, 0.01);
   nekFilter::setup("pr",  nrs->meshV, 2, 0.1);
   nekFilter::setup("ttt", cds->mesh[0], 1, 0.5);
   ```

6. Apply filter in `UDF_ExecuteStep`. Use tag to apply the specific filter.
   ```
   nekFilter::setup("vel", o_UX);
   nekFilter::setup("vel", o_UY);
   nekFilter::setup("vel", o_UZ);
   ```

### Verification

We use Nek5000's filter as a reference.

Every `nrs->isCheckpointStep=1`, we call `userchk` twice. 
At the first call, it copies the solutions to work arrays and applies filter in USR.
Then, filter is applied in UDF and the second `userchk` compare results from both.

We expect to see the folloeing in the logfile.

- NekRS
  ```
  filt trn (rs)   vel: 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 0.9000
  filt trn (rs)  temp: 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 0.9975 0.9900
  ```

- Nek5000
  ```
  # first set of filter
  filt amp 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.1000
  filt trn 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 0.9000

  # second set of filter
  filt amp 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0025 0.0100
  filt trn 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 0.9975 0.9900
  ```

- max. abs. err
  ```
  # step             err_vx      err_vy      err_vz      err_pr    err_temp      err_s1
       2 qfilt:  2.2751E-10  2.2732E-10  2.2727E-10  1.9351E-11  7.4273E-08  4.3079E-10
       4 qfilt:  2.2733E-10  2.2698E-10  2.2688E-10  1.9431E-11  1.3834E-07  1.2380E-09
       6 qfilt:  2.3557E-10  2.3501E-10  2.3485E-10  6.9847E-09  1.9508E-07  2.3552E-09
       8 qfilt:  2.3513E-10  2.3462E-10  2.3471E-10  1.6851E-09  2.4363E-07  3.6342E-09
      10 qfilt:  2.3316E-10  2.3254E-10  2.3235E-10  8.9727E-10  2.8512E-07  4.9855E-09
  ```

