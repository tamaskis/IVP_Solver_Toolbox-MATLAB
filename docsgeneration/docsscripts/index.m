%% ODE Solver Toolbox Documentation
%
% <<ode_solver_toolbox.png>>
%
% Copyright © 2021 Tamas Kis
%% Technical Documentation
% Click <https://tamaskis.github.io/documentation/Fixed-Step_ODE_Solvers.pdf here>.
%% Methology
% All of the ODE solvers in this toolbox are implemented so that they can
% be used individually without relying on any other function. Additionally,
% all of the ODE solvers support both "time detection" (solving until some
% final time) _and_ "event detection" (solving until some event occurs).
%% Installation
% The toolbox can be downloaded from <TODO: file exchange link> or 
% <TODO: github link>.
% The downloaded zip folder contains the following:
%
% * *docs* → Contains the HTML documentation. To open a copy of the HTML
% documentation locally on your computer (without need of an internet
% connection), open docs/index.html.
% * *examples* → Contains examples for using various functions of the ODE solver toolbox, as well as examples for some more elementary concepts discussed in the <https://tamaskis.github.io/documentation/Fixed-Step_ODE_Solvers.pdf technical documentation>.
% * *INSTALL* → Contains the toolbox installer (ODE Solver Toolbox.mltbx).
% * *licenses* → Contains the software licenses.
% * *README.md* → Markdown documentation for GitHub repository.
% * *Technical Documentation* → Contains the technical documentation (Fixed_Step_ODE_Solvers.pdf).
% * *toolbox* → Contains all the functions specific to this toolbox.
% * *toolbox/lib* → External libraries/functions required by this toolbox.
%
% *To install as a toolbox*, simply open "ODE Solver Toolbox.mltbx" in the 
% "INSTALL" folder. MATLAB will automatically perform the installation and 
% add all the functions included in the toolbox to the MATLAB search path.
%
% Alternatively, all the functions in the "toolbox" folder can be used
% independently, _with the exception_ of the functions in the
% "toolbox/generatingequations" folder, which require the
% |same_denominator| function in "toolbox/lib".
%% Explicit Runge-Kutta (Single-Step) Methods
% * <RK1_euler_doc.html *|RK1_euler|*> Euler method (1st-order).
% * <RK2_doc.html *|RK2|*> Midpoint method (2nd-order).
% * <RK2_heun_doc.html *|RK2_heun|*> Heun's second-order method (2nd-order).
% * <RK2_ralston_doc.html *|RK2_ralston|*> Ralston's second-order method (2nd-order).
% * <RK3_doc.html *|RK3|*> (Kutta's) Runge-Kutta third-order method (3rd-order).
% * <RK3_heun_doc.html *|RK3_heun|*> Heun's third-order method (3rd-order).
% * <RK3_ralston_doc.html *|RK3_ralston|*> Ralston's third-order method (3rd-order).
% * <SSPRK3_doc.html *|SSPRK3|*> Strong stability preserving Runge-Kutta third-order method (3rd-order).
% * <RK4_doc.html *|RK4|*> (Classic) Runge-Kutta fourth-order method (4th-order).
% * <RK4_ralston_doc.html *|RK4_ralston|*> Ralston's fourth-order method (4th-order).
% * <RK4_38_doc.html *|RK4_38|*> 3/8-rule fourth-order method (4th-order).
%% Adams-Bashforth (Multistep Predictor) Methods
% * <AB2_doc.html *|AB2|*> Adams-Bashforth 2nd-order method.
% * <AB3_doc.html *|AB3|*> Adams-Bashforth 3rd-order method.
% * <AB4_doc.html *|AB4|*> Adams-Bashforth 4th-order method.
% * <AB5_doc.html *|AB5|*> Adams-Bashforth 5th-order method.
% * <AB6_doc.html *|AB6|*> Adams-Bashforth 6th-order method.
% * <AB7_doc.html *|AB7|*> Adams-Bashforth 7th-order method.
% * <AB8_doc.html *|AB8|*> Adams-Bashforth 8th-order method.
%% Adams-Bashforth-Moulton (Multistep Predictor-Corrector) Methods
% * <ABM2_doc.html *|ABM2|*> Adams-Bashforth-Moulton 2nd-order method.
% * <ABM3_doc.html *|ABM3|*> Adams-Bashforth-Moulton 3rd-order method.
% * <ABM4_doc.html *|ABM4|*> Adams-Bashforth-Moulton 4th-order method.
% * <ABM5_doc.html *|ABM5|*> Adams-Bashforth-Moulton 5th-order method.
% * <ABM6_doc.html *|ABM6|*> Adams-Bashforth-Moulton 6th-order method.
% * <ABM7_doc.html *|ABM7|*> Adams-Bashforth-Moulton 7th-order method.
% * <ABM8_doc.html *|ABM8|*> Adams-Bashforth-Moulton 8th-order method.
%% One-Step Propagation
% * <RK1_euler_step_doc.html *|RK1_euler_step|*> Propagates the state vector forward one time step using the Euler method (1st-order).
% * <RK2_step_doc.html *|RK2_step|*> Propagates the state vector forward one time step using the midpoint method (2nd-order).
% * <RK2_heun_step_doc.html *|RK2_heun_step|*> Propagates the state vector forward one time step using Heun's second-order method (2nd-order).
% * <RK2_ralston_step_doc.html *|RK2_ralston_step|*> Propagates the state vector forward one time step using Ralston's second-order method (2nd-order).
% * <RK3_step_doc.html *|RK3_step|*> Propagates the state vector forward one time step using (Kutta's) Runge-Kutta third-order method (3rd-order).
% * <RK3_heun_step_doc.html *|RK3_heun_step|*> Propagates the state vector forward one time step using Heun's third-order method (3rd-order).
% * <RK3_ralston_step_doc.html *|RK3_ralston_step|*> Propagates the state vector forward one time step using Ralston's third-order method (3rd-order).
% * <SSPRK3_step_doc.html *|SSPRK3_step|*> Propagates the state vector forward one time step using the strong stability preserving Runge-Kutta third-order method (3rd-order).
% * <RK4_step_doc.html *|RK4_step|*> Propagates the state vector forward one time step using the (classic) Runge-Kutta fourth-order method (4th-order).
% * <RK4_ralston_step_doc.html *|RK4_ralston_step|*> Propagates the state vector forward one time step using Ralston's fourth-order method (4th-order).
% * <RK4_38_step_doc.html *|RK4_38_step|*> Propagates the state vector forward one time step using the 3/8-rule Runge-Kutta fourth-order method (4th-order).
%% Tools for Matrix-Valued ODEs
% * <Matrix_ODE_Example_doc.html *Matrix ODE Example*> Example for solving a matrix-valued ODE.
% * <odefun_mat2vec_doc.html *|odefun_mat2vec|*> Transforms a matrix-valued ODE into a vector-valued ODE.
% * <odeIC_mat2vec_doc.html *|odeIC_mat2vec|*> Transforms the initial condition for a matrix-valued ODE into the initial condition for the corresponding vector-valued ODE.
% * <odesol_vec2mat_doc.html *|odesol_vec2mat|*> Transforms the solution matrix for a vector-valued ODE into the solution array for the corresponding matrix-valued ODE.
%% Generating ODE Solver Equations
% * <AB_coefficients_doc.html *|AB_coefficients|*> Coefficients for the mth-order Adams-Bashforth predictor.
% * <AM_coefficients_doc.html *|AM_coefficients|*> Coefficients for the mth-order Adams-Moulton corrector.
% * <AB_predictor_doc.html *|AB_predictor|*> mth-order Adams-Bashforth predictor.
% * <AM_corrector_doc.html *|AM_corrector|*> mth-order Adams-Moulton corrector.
% * <ABM_equations_doc.html *|ABM_equations|*> mth-order Adams-Bashforth-Moulton equations.
% * <tableau2eqns_doc.html *|tableau2eqns|*> Propagation equations from Butcher tableau for explicit Runge-Kutta methods.
%% External Libraries
% * <https://www.mathworks.com/matlabcentral/fileexchange/95618-convert-fractions-to-same-denominator-same_denominator *|same_denominator|*> Scales a set of fractions so they each have the same denominator.