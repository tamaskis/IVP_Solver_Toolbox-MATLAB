# IVP Solver Toolbox [![View ODE Solver Toolbox on File Exchange](https://www.mathworks.com/matlabcentral/images/matlab-file-exchange.svg)](https://www.mathworks.com/matlabcentral/fileexchange/103975-ode-solver-toolbox)

Collection of fixed-step IVP solvers, functions to generate IVP solver equations, and functions for transforming matrix IVPs to vector IVPs.


## Documentation

[Toolbox Documentation](https://tamaskis.github.io/IVP_Solver_Toolbox-MATLAB/)\
[Technical Documentation](https://tamaskis.github.io/documentation/Solving_Initial_Value_Problems_for_ODEs.pdf)

To open the home page of the toolbox documentation in MATLAB, type

```
doc_IST
```

in the Command Window. To open the documentation of a specific function with name `function_name` from the Command Window, type

```
doc_IST function_name
```

To open the PDF file with the technical documentation (Solving_Initial_Value_Problems_for_ODEs.pdf) from the Command Window, type

```
doc_IST tech
```

## Main IVP Solver Function



## Explicit Runge-Kutta (Single-Step) Methods

`RK1_euler`\
`RK2`\
`RK2_heun`\
`RK2_ralston`\
`RK3`\
`RK3_heun`\
`RK3_ralston`\
`SSPRK3`\
`RK4`\
`RK4_ralston`\
`RK4_38`


## Adams-Bashforth (Multistep Predictor) Methods
`AB2`\
`AB3`\
`AB4`\
`AB5`\
`AB6`\
`AB7`\
`AB8`


## Adams-Bashforth-Moulton (Multistep Predictor-Corrector) Methods
`ABM2`\
`ABM3`\
`ABM4`\
`ABM5`\
`ABM6`\
`ABM7`\
`ABM8`


## One-Step Propagation

`RK1_euler_step`\
`RK2_step`\
`RK2_heun_step`\
`RK2_ralston_step`\
`RK3_step`\
`RK3_heun_step`\
`RK3_ralston_step`\
`SSPRK3_step`\
`RK4_step`\
`RK4_ralston_step`


## Tools for Matrix-Valued IVPs
`odefun_mat2vec`\
`ivpIC_mat2vec`\
`ivpsol_vec2mat`


## Generating IVP Solver Equations
`AB_coefficients`\
`AM_coefficients`\
`AB_predictor`\
`AM_corrector`\
`ABM_equations`\
`tableau2eqns`
