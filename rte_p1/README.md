# (WIP) Radiation P1 Model

- version: v24 (repo/next, 08/06/24)

### Equations

- Heat transfer equation   
  $$\frac{\partial T}{\partial t} + \vec{u}\cdot \nabla T = \nabla\cdot(\frac{1}{RePr}\nabla T - \frac{4\tau}{B_o}(T^4 - G)$$   
  with Neumann BC: $q\cdot\vec{n} = -\frac{4}{3\tau}{B_o}\nabla G\cdot\vec{n} = -\frac{2\epsilon}{B_o (2-\epsilon)}(G-T^4)$

- Incident Radiation intensity   
  $$-\nabla\cdot\frac{1}{3\tau}\nabla G + T^4$$   
  with Robin BC: $\frac{1}{3\tau^2}(\nabla G\cdot\vec{n}) = \frac{\epsilon}{2(2-\epsilon)\tau} (G-T^4)$


