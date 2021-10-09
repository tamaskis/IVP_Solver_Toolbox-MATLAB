%% ODE Solver Toolbox Documentation
%
% <<ode_solver_toolbox.png>>
%
% Copyright © 2021 Tamas Kis
%% EXAMPLES
% Since all of the functions have an identical syntax, only one set of
% examples is included. These examples use the <RK4_doc.html |RK4|> function, 
% but could also be solved using any of the other solver functions in this 
% toolbox.
%% Single-Step Methods
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
%% Multistep Predictor Methods
%% Multistep Predictor-Corrector Methods
% * <ABM8_doc.html *|ABM8|*> Adams-Bashforth-Moulton 8th-order method.
%% Generating ODE Solver Equations
% * <AB_coefficients_doc.html *|AB_coefficients|*> Coefficients for the nth-order Adams-Bashforth predictor.
% * <AM_coefficients_doc.html *|AM_coefficients|*> Coefficients for the nth-order Adams-Moulton corrector.
% * <AB_predictor_doc.html *|AB_predictor|*> nth-order Adams-Bashforth predictor.
% * <AM_corrector_doc.html *|AM_corrector|*> nth-order Adams-Moulton corrector.
% * <ABM_equations_doc.html *|ABM_equations|*> nth-order Adams-Bashforth-Moulton equations.
% * <tableau2eqns_doc.html *|tableau2eqns|*> Propagation equations from Butcher tableau for explicit Runge-Kutta methods.