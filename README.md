This package consists of functionalities useful for basic engineering studies.

Current Implementation : Finite Element Analysis on 2D problems has been implemented.\
In Progress : Literally everything.

[Github Repo Link](https://github.com/sayanbiswas023/mechy)\
[Pypi link](https://pypi.org/project/mechy/)

#INSTALLATION GUIDE:

```
pip install mech
```
```
pip3 install mechy
```

#DOCUMENTATION:

##2D Finite Element Analysis
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


```
plot_type='all'
```
Variable plot_type should be among the list of possible output demands :
``` ['u1','u2','s11','s22','s12','e11','e22','e12','all'] 
```

The input file format is similar to standard simulating softwares as abaqus

```
### Input file syntax:

    *Material
    500, 0.06             # E and nu
    
    *Node
    1, 0.0, 0.0
    2, 0.0, 1.0
    3, 1.0, 1.0
    4, 1.0, 0.0
    
    *Element
    1, 1, 2, 3, 4        # elemId, nodeid of vertices
    
    *Step
    *Boundary
    1, 1, 2, 0.0          # nodeId, dof1, dof2, value
    2, 1, 1, 0.0
    3, 1, 1, 0.01
    4, 1, 1, 0.01
    4, 2, 2, 0.0
```

<br>

<img src="https://media.tenor.com/rJ4wiVqdmcAAAAAC/sponge-bob-hammer.gif" width="500" height="400">
