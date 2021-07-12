%% |tableau2eqns|
% Propagation equations from Butcher tableau.
% 
% <ODE_Solver_Toolbox_Contents.html Back to ODE Solver Toolbox Contents>.
%% Syntax
%   tableau2eqns(T)
%   tableau2eqns(T,name)
%   tableau2eqns(__,'decimal')
%% Input/Output Parameters
% <html>
%   <table border=1>
%       <tr>
%           <td></td>
%           <td style="text-align:center"><b>Symbol</b></td>
%           <td style="text-align:center"><b>Variable</b></td>
%           <td style="text-align:center"><b>Description</b></td>
%           <td style="text-align:center"><b>Format</b></td>
%       </tr>
%       <tr>
%           <td rowspan="3" style="text-align:center"><b>Input</b></td>
%           <td style="text-align:center">-</td>
%           <td><TT>T</TT></td>
%           <td>Butcher tableau</td>
%           <td style="text-align:center">1×FIX SIZE<BR>double</td>
%       </tr>
%       <tr>
%           <td style="text-align:center">-</td>
%           <td><TT>name</TT></td>
%           <td>(OPTIONAL) name of Runge-Kutta method corresponding to the Butcher tableau</td>
%           <td style="text-align:center">string</td>
%       </tr>
%       <tr>
%           <td style="text-align:center">-</td>
%           <td><TT>type</TT></td>
%           <td>(OPTIONAL) print coefficients as <TT>'decimal'</TT> or <TT>'fraction'</TT> (defaults to <TT>'fraction'</TT>)</td>
%           <td style="text-align:center">char</td>
%       </tr>
%   </table>
% </html>
%% Example #1: RK2 (Midpoint method)
T = [0    0    0;
     1/2  1/2  0;
     0    0    1];
tableau2eqns(T,"RK2 (Midpoint method)");
%% Example #2: RK2 (Heun's method)
T = [0  0    0;
     1  1    0;
     0  1/2  1/2];
tableau2eqns(T,"RK2 (Heun's method)");
%% Example #3: (Ralston's method)
T = [0    0    0;
     2/3  2/3  0;
     0    1/4  3/4];
tableau2eqns(T,"RK2 (Ralston's method)");
%% Example #4: (Kutta's third-order method)
T = [0     0    0    0;
     1/2   1/2  0    0;
     1    -1    2    0;
     0     1/6  2/3  1/6];
tableau2eqns(T,"RK3 (Kutta's third-order method)");
%% Example #5: (Heun's third-order method)
T = [0    0    0    0;
     1/3  1/3  0    0;
     2/3  0    2/3  0;
     0    1/4  0    3/4];
tableau2eqns(T,"RK3 (Heun's third-order method)");
%% Example #6: (Ralston's third-order method)
T = [0    0    0    0;
     1/2  1/2  0    0;
     3/4  0    3/4  0;
     0    2/9  1/3  4/9];
tableau2eqns(T,"RK3 (Ralston's third-order method)");
%% Example #7: (Strong Stability Preserving Runge-Kutta third-order method)
T = [0    0     0    0;
     1    1     0    0;
     1/2  1/4   1/4  0;
     0    1/6   1/6  2/3];
tableau2eqns(T,"SSPRK3 (Strong Stability Preserving Runge-Kutta third-order method)");
%% Example #8: (Runge-Kutta fourth-order method)
T = [0     0     0    0    0;
     1/2   1/2   0    0    0;
     1/2   0     1/2  0    0;
     1     0     0    1    0;
     0     1/6   1/3  1/3  1/6];
tableau2eqns(T,"RK4 (Runge-Kutta fourth-order method)");
%% Example #9: (Ralston's fourth-order method)
T = [0           0            0             0            0;
     0.4         0.4          0             0            0;
     0.45573725  0.29697761   0.15875964    0            0;
     1           0.21810040  -3.05096516    3.83286476   0;
     0           0.17476028  -0.55148066    1.20553560   0.17118478];
tableau2eqns(T,"RK4 (Ralston's fourth-order method)",'decimal');
%% Example #10: (3/8-Rule fourth-order method)
T = [0     0     0    0    0;
     1/3   1/3   0    0    0;
     2/3  -1/3   1    0    0;
     1     1    -1    1    0;
     0     1/8   3/8  3/8  1/8];
tableau2eqns(T,"3/8-Rule fourth-order method");