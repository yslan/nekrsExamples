# Del operators

Verifying and test derivaives operators in `opSEM::`

- Version: v24pre, (repo/next cloned at Oct. 19, 2024)

- Domain: 3D Flow past a sphere inside cylinder container

  with dimension:

   - sphere radius = 1

   - cylinder radius = 4.4

  The mesh comes from my [OneSphere](https://github.com/yslan/OneSphere) project.
 
## Tests:

1. volume

2. surface area

3. divergence

   The current NekRS version has a bug in `auto o_tdiv = opSEM::strongDivergence`.
   One will get infinite loop.A
   Use another interface: `void opSEM::strongDivergence(... , o_out);` instead.

   - Input:   
     $$\vec{f} = [x^2\sin(y), y\cos(z), \exp(x)\sin(z)]$$

   - Output:   
     $$\nabla\cdot\vec{f}$$

   - Exact solution:    
     $$2x\sin(y) + \cos(z) + \exp(x)\cos(z)$$

4. gradient

   - Input:    
     $$f(x,y,z) = \exp(-x^2-y^2) * \sin(z)$$

   - Output:   
     $$\nabla f$$

   - Exact solution: 
     $$[-2x\sin(z)\exp(-x^2-y^2), -2y\sin(z)\exp(-x^2-y^2), \cos(z)\exp(-x^2-y^2)]$$

5. laplacian

   Due to the bug in 4 and `o_out = strongDivergence` will cause issue.
   Here, use the inline function `strongLaplacian` as a workaround.

   - Input:    
     $$f(x,y,z) = \exp(-x^2-y^2) * \sin(z)$$

   - Output:   
     $$\nabla\cdot\left(a \nabla f\right)$$
     where $a(x,y,z) = 1+z$
     Note that, there is no minus sign at front.

   - Exact solution:    
     $$\exp(-x^2-y^2)\left[(1+z)(4x^2+4y^2-5)\sin(z) + \cos(z)\right]$$

(WIP) Expected output:
```
div chk N=  9 err= 6.47721125E-08 mag= 1.71421498E+01 4.40000000E+00 8.14412143E+01(5k)
div chk N=  9 err= 6.67595580E-08 mag= 1.71421498E+01 4.40000000E+00 8.14412143E+01(5k)
vol chk N=  9 err=  9.95243679E-14 (5k)
area chk N=  9 errs=  2.82715972E-16 2.74148821E-16 5.84123908E-16 3.50474345E-16 (5k)
vol chk N= 9 err= 4.88232748e-15
div chk N= 9 err (inf,L2) = 3.72665209e-07 6.47721154e-08
grad chk N= 9 err (inf,L2) = 1.57675467e-08 1.94239270e-09
lap chk N= 9 err (inf,L2) = 3.89993326e-07 6.43594901e-08
lap chk N= 9 err (inf,L2) = 4.45114789e-07 9.33265581e-08
```
