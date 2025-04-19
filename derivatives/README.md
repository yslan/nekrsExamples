# opSEM playground

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
     $$\nabla\cdot\left(a \nabla f)$$
     where $a(x,y,z) = 1+z$
     Note that, there is no minus sign at front.

   - Exact solution:    
     $$\exp(-x^2-y^2)\left[(1+z)(4x^2+4y^2-5)\sin(z) + \cos(z)\right]$$

