# Example of Robin BC

This example show how to solve scalars with Robin type BC, to be more specific, the Newton cooling BC `"c  "` in Nek5000.

- version: v24 (repo/next 07/10/24 + some fixes that will be released later)

### Implementation

Use Neumann BC `"f  "`, tag boundary ID for Robin in `UDF_Setup`. 
The RHS (derivative part) can be computed in `.oudf` and the LHS will be added into `h2` via the implicit linear term

- Newton cooling   
  $$- \kappa \nabla T \cdot {\vec n} = h_c (T - T_{\infty})$$

- Robin
  $$\alpha \nabla T \cdot {\vec n} + \beta T = \gamma$$      
  Therefore, $\gamma = h_c T_{\infty}$, $\beta = h_c$, $\alpha=\kappa$


### Verification

See the [Nek5000 example](https://github.com/Nek5000/NekExamples/tree/master/robin). There is a PDF for the derivation.     
The averaged temperature has theoretical exponential decay rate `exp(-lamb*time)`.     
In `UDF_ExecuteStep`, we print tbar and the `ratio = theo. rate / tbar` will conv. to a  constant.

```
$ grep lamb logfile |tail


 #     step   time         tbar              lambda       theo. rate   rate / tbar -> const.
 tbar: 1991   9.9550e+00   9.3359e-06 lamb   1.1597e+00   9.6902e-06   1.0380e+00
 tbar: 1992   9.9600e+00   9.2819e-06 lamb   1.1597e+00   9.6342e-06   1.0380e+00
 tbar: 1993   9.9650e+00   9.2283e-06 lamb   1.1597e+00   9.5785e-06   1.0380e+00
 tbar: 1994   9.9700e+00   9.1749e-06 lamb   1.1597e+00   9.5231e-06   1.0380e+00
 tbar: 1995   9.9750e+00   9.1219e-06 lamb   1.1597e+00   9.4681e-06   1.0380e+00
 tbar: 1996   9.9800e+00   9.0691e-06 lamb   1.1597e+00   9.4133e-06   1.0380e+00
 tbar: 1997   9.9850e+00   9.0167e-06 lamb   1.1597e+00   9.3589e-06   1.0380e+00
 tbar: 1998   9.9900e+00   8.9646e-06 lamb   1.1597e+00   9.3048e-06   1.0380e+00
 tbar: 1999   9.9950e+00   8.9127e-06 lamb   1.1597e+00   9.2510e-06   1.0380e+00
 tbar: 2000   1.0000e+01   8.8612e-06 lamb   1.1597e+00   9.1975e-06   1.0380e+00
```
