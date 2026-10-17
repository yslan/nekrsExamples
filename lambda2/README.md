## Lambda2

Compute lambda2 on device.

Tested with NekRS v26


### Setup

Copy the hpp to the case filder and include in udf both inside and outside the
`__okl__` block. It follows `nrs->Qcriterion`. See `lmabda2.hpp` for more interfaces.

```
#include "lambda2.hpp"
#ifdef __okl__ 
#include "lambda2.hpp" // kernel
#endif
```

```
auto o_lambda2 = lambda2(false); // F = no gsavg
```


### TODO
- fp64 compareison between Nek5000 one
- explicit filter
- merge into NekRS
