# IVP Solver Toolbox [![View IVP Solver Toolbox on File Exchange](https://www.mathworks.com/matlabcentral/images/matlab-file-exchange.svg)](https://www.mathworks.com/matlabcentral/fileexchange/103975-ivp-solver-toolbox)

Collection of fixed-step IVP solvers, functions to generate IVP solver equations, and functions for transforming matrix IVPs to vector IVPs.


## Documentation

[Toolbox Documentation](https://tamaskis.github.io/IVP_Solver_Toolbox-MATLAB/)\
[Technical Documentation](https://tamaskis.github.io/files/Solving_Initial_Value_Problems_for_ODEs.pdf)

To open the home page of the toolbox documentation in MATLAB, type

```
doc_IVP
```

in the Command Window. To open the documentation of a specific function with name `function_name` from the Command Window, type

```
doc_IVP function_name
```

To open the PDF file with the technical documentation (Solving_Initial_Value_Problems_for_ODEs.pdf) from the Command Window, type

```
doc_IVP tech
```

## Main IVP Solver Function

`[t,y] = solve_ivp(f,[t0,tf],y0,h)`\
`[t,y] = solve_ivp(f,{t0,C},y0,h)`\
`[t,y] = solve_ivp(__,method)`\
`[t,y] = solve_ivp(__,method,wb)`


## Matrix-Valued IVP Solver Function

`[t,M] = solve_ivp_matrix(F,[t0,tf],M0,h)`\
`[t,M] = solve_ivp_matrix(F,{t0,C},M0,h)`\
`[t,M] = solve_ivp_matrix(__,p,method)`\
`[t,M] = solve_ivp_matrix(__,p,method,wb)`


## Utilities for IVP Solvers
`expand_ivp_arrays`\
`mat2vec_fun`\
`mat2vec_IC`\
`mat2vec_C`\
`vec2mat_sol`


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


## Generating IVP Solver Equations
`AB_coefficients`\
`AM_coefficients`\
`AB_predictor`\
`AM_corrector`\
`ABM_equations`\
`tableau2eqns`