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
% * *INSTALL* → Contains the toolbox installer (ODE Solver Toolbox.mltbx).
% * *licenses* → Contains the software licenses.
% * *README.md* → Markdown documentation for GitHub repository.
% * *Technical Documentation* → Contains the technical documentation (Fixed_Step_ODE_Solvers.pdf).
% * *tests* → Code for unit tests.
% * *toolbox* → Contains all the functions specific to this toolbox.
% * *toolbox/lib* → External libraries/functions required by this toolbox.
%
% *To install as a toolbox*, simply open "ODE Solver Toolbox.mltbx" in the 
% "INSTALL" folder. MATLAB will automatically perform the installation and 
% add all the functions included in the toolbox to the MATLAB search path.
%
% Alternatively, all the functions in the "toolbox/rungekutta", 
% "toolbox/adamsbashforthmoulton", and "toolbox/adamsbashforth" folders can
% be used independently.
%% Explicit Runge-Kutta (Single-Step) Methods
% * <euler_doc.html *|euler|*> Euler method (1st-order).
% * <RK2_doc.html *|RK2|*> Midpoint method (2nd-order).
% * <RK2_heun_doc.html *|RK2_heun|*> Heun's second-order method (2nd-order).
% * <RK2_ralston_doc.html *|RK2_ralston|*> Ralston's second-order method (2nd-order).
% * <RK3_doc.html *|RK3|*> (Kutta's) Runge-Kutta third-order method (3rd-order).
% * <RK3_heun_doc.html *|RK3_heun|*> Heun's third-order method (3rd-order).
% * <RK3_ralston_doc.html *|RK3_ralston|*> Ralston's third-order method (3rd-order).
% * <SSPRK3_doc.html *|SSPRK3|*> Strong stability preserving Runge-Kutta third-order method (3rd-order).
% * <RK4_doc.html *|RK4|*> (Classic) Runge-Kutta fourth-order method (4th-order).
% * <RK4_ralston_doc.html *|RK4_ralston|*> Ralston's fourth-order method (4th-order).
% * <RK4_38_doc.html *|RK4_38|*> 3/8-Rule fourth-order method (4th-order).
%% Adams-Bashforth (Multistep Predictor) Methods
%% Adams-Bashforth-Moulton (Multistep Predictor-Corrector) Methods
% * <ABM2_doc.html *|ABM2|*> Adams-Bashforth-Moulton 2nd-order method.
% * <ABM3_doc.html *|ABM3|*> Adams-Bashforth-Moulton 3rd-order method.
% * <ABM4_doc.html *|ABM4|*> Adams-Bashforth-Moulton 4th-order method.
% * <ABM5_doc.html *|ABM5|*> Adams-Bashforth-Moulton 5th-order method.
% * <ABM6_doc.html *|ABM6|*> Adams-Bashforth-Moulton 6th-order method.
% * <ABM7_doc.html *|ABM7|*> Adams-Bashforth-Moulton 7th-order method.
% * <ABM8_doc.html *|ABM8|*> Adams-Bashforth-Moulton 8th-order method.
%% Generating ODE Solver Equations
% * <AB_coefficients_doc.html *|AB_coefficients|*> Coefficients for the mth-order Adams-Bashforth predictor.
% * <AM_coefficients_doc.html *|AM_coefficients|*> Coefficients for the mth-order Adams-Moulton corrector.
% * <AB_predictor_doc.html *|AB_predictor|*> mth-order Adams-Bashforth predictor.
% * <AM_corrector_doc.html *|AM_corrector|*> mth-order Adams-Moulton corrector.
% * <ABM_equations_doc.html *|ABM_equations|*> mth-order Adams-Bashforth-Moulton equations.
% * <tableau2eqns_doc.html *|tableau2eqns|*> Propagation equations from Butcher tableau for explicit Runge-Kutta methods.
%% External Libraries
% * <https://www.mathworks.com/matlabcentral/fileexchange/95618-convert-fractions-to-same-denominator-same_denominator *|same_denominator|*> Scales a set of fractions so they each have the same denominator.