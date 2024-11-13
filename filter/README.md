# Explicit filter

Here we apply the "explicit" filter used in Nek5000.
The case is built on the top of the either example.
TODO: maybe change to an example that does need filtering.

- Version: NekRS v24-pre (repo/next, tested on 08/18/24)
- `v23/` is added.


### Formulation

See the Nek5000 [docuentation](https://nek5000.github.io/NekDoc/problem_setup/filter.html).

Basically, in each direction, it changes the basis, rescale the magnitude with diagonal matrix then change the basis back. The whole action is stored into a 1D matrix and the kernel simply do max-vec in each directions, elements by elements.


### Case setup

1. Put the files `nekFilter.hpp` and `nekFilter.okl` to the case folder.

2. `#include "nekFilter.hpp"` in UDF file.

3. Turn off `regularization` for the desired field.       
   This is user's responsibility.
   We don't check this in the code since user may want to apply filter at will to arbitrary fields, says post-processing the Q-criterion.

4. Call `nekFilter::buildKernel(kernelInfo);` in the UDF function `UDF_LoadKernels`.

5. Call the setup function with a name tag of the filter, following by the mesh, filterNode and filterWeight. It can have multiple filters for each fields. e.g.
   ```
   nekFilter::setup("vel", nrs->meshV, 1, 0.01);
   nekFilter::setup("pr",  nrs->meshV, 2, 0.1);
   nekFilter::setup("ttt", cds->mesh[0], 1, 0.5);
   ```

6. Apply filter in `UDF_ExecuteStep`. Use the name tag to apply the specific filter.
   ```
   nekFilter::setup("vel", o_UX);
   nekFilter::setup("vel", o_UY);
   nekFilter::setup("vel", o_UZ);

   nekFilter::setup("ttt", cds->o_S);
   ```

### Verification

We use Nek5000's filter as a reference.

Every `nrs->isCheckpointStep=1`, we call `userchk` twice. 
At the first call, it copies the solutions to work arrays and applies filter in Nek5000.
Then, NekRS' filter is applied after and the second `userchk` compares results from both.

We expect to see the folloeing in the logfile.

- NekRS
  ```
  filt trn (rs)   vel: 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 0.9889 0.9556 0.9000
  filt trn (rs)  temp: 1.0000 1.0000 1.0000 1.0000 1.0000 0.9996 0.9984 0.9964 0.9936 0.9900
  ```

- Nek5000
  ```
  # first set of filter
  filt amp 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0111 0.0444 0.1000
  filt trn 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 0.9889 0.9556 0.9000

  # second set of filter
  filt amp 0.0000 0.0000 0.0000 0.0000 0.0000 0.0004 0.0016 0.0036 0.0064 0.0100
  filt trn 1.0000 1.0000 1.0000 1.0000 1.0000 0.9996 0.9984 0.9964 0.9936 0.9900
  ```

- Relative Linf err
  ```
  # step             err_vx      err_vy      err_vz      err_pr    err_temp      err_s1
       2 qfilt:  1.2251E-15  9.6470E-16  1.0243E-15  1.2153E-15  1.2251E-15  1.2252E-15
       4 qfilt:  1.2249E-15  9.6308E-16  1.0232E-15  1.3892E-15  1.4291E-15  1.4292E-15
       6 qfilt:  1.4288E-15  1.1538E-15  1.2492E-15  1.1121E-15  1.2247E-15  1.2248E-15
       8 qfilt:  1.2244E-15  9.5988E-16  1.0210E-15  1.1101E-15  1.2244E-15  1.2246E-15
      10 qfilt:  1.2242E-15  9.5829E-16  1.0199E-15  1.1080E-15  1.2242E-15  1.2244E-15
  ```


### Single file

See the case under `single_file/` with all functions put into the same file.

See also `single_zero_file/` for removing the ethier setup and verification.



