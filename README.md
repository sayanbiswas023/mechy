Here is a demo code to display all outputs of a 2D Finite Element Analysis on a dogbone.

```
from mechy import FEM
FILE_PATH='test.txt'
plot_type='all'

s=FEM()
s.fem2d(FILE_PATH,plot_type)
```
##Displacement_X
![Displacement_X](./mechy/images/u1.png)
##Stress_11
![stress_11](./mechy/images/s11.png)
##stress_22
![stress_22](./mechy/images/s22.png)
##strain_12
![strain_12](./mechy/images/e12.png)
