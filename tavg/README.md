# Modifed tavg

- Version, NekRS v24pre (repo/next, cloned at Oct. 19)

- Setup:
  Replace the following two files:
  ```
  nekrs/src/plugins/tavg.cpp
  nekrs/src/plugins/tavg.hpp
  ```

## Usage and Old Features

The NekRS's default still works (See turbPipe example).
On the top of that, two features are added.


1. Add auto setup for `avg_all`
   Previously, user has to manually set something like this
   ```
   tavgFields.push_back({o_u});
   tavgFields.push_back({o_v});
   tavgFields.push_back({o_w});
   tavgFields.push_back({o_temp});

   tavg::setup(nrs->fieldOffset, tavgFields);
   ```
   to have `E[X]`. And, tavg will save every fields into scalars slot of the file `tavg0.f0000*`

   Now, we add a new mode:
   ```
   tavg::setup(nrs);
   ```

   It automatically does the same thing as `avg_all`.
   ```
   // E[X] into "avg"
   tavgFields.push_back({o_u});
   tavgFields.push_back({o_v});
   tavgFields.push_back({o_w});
   tavgFields.push_back({o_p});
   tavgFields.push_back({o_temp}); // if present
   tavgFields.push_back({o_S1});   // if present

   // E[X^2] into "rms"
   tavgFields.push_back({o_u, o_u});
   tavgFields.push_back({o_v, o_v});
   tavgFields.push_back({o_w, o_w});
   tavgFields.push_back({o_p, o_p});
   tavgFields.push_back({o_temp, o_temp}); // if present
   tavgFields.push_back({o_S1, o_S1});     // if present

   // E[XY] into "rm2"
   tavgFields.push_back({o_u, o_v});
   tavgFields.push_back({o_v, o_w});
   tavgFields.push_back({o_w, o_u});
   ```

2. Add "FP64" and "reset" to `tavg::outfld(mesh, FP64 = true, reset = true)`.
   The default values for both modes are `true` which dumps files in double precision and reset the `atime` everytime it dumpes a file.
   Now, it supports
   ```
   // default
   tavg::outfld(mesh);

   // FP64, no-reset
   tavg::outfld(mesh, true, false);

   // FP32, no-reset
   tavg::outfld(mesh, false, false);
   ```

   The reason for FP64 is, Paraview's new Nek5000 reader only supports single prevision. 
   https://gitlab.kitware.com/paraview/paraview/-/issues/22747
