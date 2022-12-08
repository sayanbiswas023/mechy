This package consists of functionalities useful for basic engineering studies.

Current Implementation : Finite Element Analysis on 2D problems has been implemented.\
In Progress : FEM problems stated in 3D.\

[Github Repo Link](https://github.com/sayanbiswas023/mechy)
INSTALLATION GUIDE:

```
pip install mech
```
```
pip3 install mechy
```

Here is a demo code to display all outputs of a 2D Finite Element Analysis on a dogbone.

```
from mechy import FEM
FILE_PATH='demo_examples/2d_test.txt'
plot_type='all'

s=FEM()
s.fem2d(FILE_PATH,plot_type)
```

![Displacement_X](https://raw.githubusercontent.com/sayanbiswas023/mechy/master/mechy/images/u1.png)
![stress_11](https://raw.githubusercontent.com/sayanbiswas023/mechy/master/mechy/images/s11.png)
![stress_22](https://raw.githubusercontent.com/sayanbiswas023/mechy/master/mechy/images/s22.png)
![strain_12](https://raw.githubusercontent.com/sayanbiswas023/mechy/master/mechy/images/e12.png)

DOCUMENTATION:
```
plot_type='all'
```
Variable plot_type should be among the list of possible output demands : ['u1','u2','s11','s22','s12','e11','e22','e12','all']
