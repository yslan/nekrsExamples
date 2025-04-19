# NekrsExamples
My own collactions of NekRS examples.     
Here we focus on porting Nek5000 functionality into NekRS

- turbOutflow        
  use nek5000's `turb_outflow` in nekrs   

- turbInflow        
  recycle velocity using findpt     

- airfoil      
  compute drag coefficient, match nek5000's `torque_calc`

- `kershaw_bps5`        
  CEED's bps5

- robin     
  Robin BC, or Newton cooling BC in Nek5000

- hmhSolver    
  solve extra Helmholtz solvers in `UDF_ExecuteStep`

- nusselt      
  Compute nusselt number on given surface

- rte_p1    
  radiation transfer, P1 model

- filter    
  Nek5000's explicit filter

- `hemi_elementFilter`
  detect CFL, only dump elelments aruond it using `elementFilter`

- `iofldFactory`
  A example of having extra iofldFactory.
  
- `tavg`
  Extend support of tavg plugin to have Nek5000's `avg_all`

- `del`
  Verify derivative operators in `opSEM::`
