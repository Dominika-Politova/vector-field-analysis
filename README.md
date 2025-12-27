## vector-field-analysis
C++ library for 3D vector field analysis with 1D Helmholtz decomposition

## Features
- 3D vector field storage on regular Cartesian grids
- 3D divergence calculation (∇·F) using central differences
- 3D curl calculation (∇×F) using central differences
- 1D Helmholtz decomposition analysis (using Helmholz Decomposition Theorem)
- 1D Poisson equation solver for potential fields (using Trapezoidal rule for integration)

## 1D Analysis Note:
The Helmholtz decomposition currently works on 1D slices (along X-axis at Y,Z center).
This provides insight into field behavior but is **not a full 3D decomposition**.

**Important limitation**: Y and Z components are placed entirely into the solenoidal part in the 1D analysis, though they may contain both potential and solenoidal contributions.

## Future Work:
- Full 3D Helmholtz decomposition implementation
- Boundary condition support
- Parallel computation support
- Hydrodynamic/Magnetic field applications

## Example of output (running example_4)
Analyzed field: F = {y, x * z, y} (100% solenoidal)

```
Divergence at (10,10,10): 0
Rotation at (10,10,10): (0.9, 0, -0.9)

    Field Analysis
 (Only x-coordinate!)

change of the rotor at several points
  i=4: rot = (0.96, 0, -0.5)
  i=5: rot = (0.95, 0, -0.5)
  i=6: rot = (0.94, 0, -0.5)
div(F_sol) in center:     0 (must be 0)

results in center (i=50, x=0.5):
  full field: (0.5, 0.25, 0.5)
  potential part: (-0, 0, 0)
  solenoidal part: (0.5, 0.25, 0.5)
  sum of the solenoidal and potential parts: (0.5, 0.25, 0.5)
  error: (0, 0, 0)
   div(F) = 0
   rot(F) = (0.5, 0, -0.5)
```
