%% |tableau2eqns|
% Propagation equations from Butcher tableau.
% 
% <index.html Back to ODE Solver Toolbox Contents>.
%% Syntax
%   tableau2eqns(T)
%   tableau2eqns(T,name)
%   tableau2eqns(__,'decimal')
%% Description
% |tableau2eqns(T)| prints the propagation equations corresponding to a
% Butcher Tableau defined by |T| to the Command Window.
%%
% |tableau2eqns(T,name)| does the same as the syntax above, but labels the
% equations with the name specified by |name|.
%%
% |tableau2eqns(__,'decimal')| does the same as the syntaxes above, but
% writes the coefficients of the propagation equations in decimal form.
%% Input/Output Parameters
% <html>
%   <table border=1>
%       <tr>
%           <td></td>
%           <td style="text-align:center"><b>Variable</b></td>
%           <td style="text-align:center"><b>Symbol</b></td>
%           <td style="text-align:center"><b>Description</b></td>
%           <td style="text-align:center"><b>Format</b></td>
%       </tr>
%       <tr>
%           <td rowspan="3" style="text-align:center"><b>Input</b></td>
%           <td style="text-align:center"><TT>T</TT></td>
%           <td style="text-align:center">-</td>
%           <td>Butcher tableau</td>
%           <td style="text-align:center">(s+1)×(s+1)<BR>double</td>
%       </tr>
%       <tr>
%           <td style="text-align:center"><TT>name</TT></td>
%           <td style="text-align:center">-</td>
%           <td>(OPTIONAL) name of Runge-Kutta method corresponding to the Butcher tableau</td>
%           <td style="text-align:center">1×1<BR>string</td>
%       </tr>
%       <tr>
%           <td style="text-align:center"><TT>type</TT></td>
%           <td style="text-align:center">-</td>
%           <td>(OPTIONAL) print coefficients as <TT>'decimal'</TT> or <TT>'fraction'</TT> (defaults to <TT>'fraction'</TT>)</td>
%           <td style="text-align:center">char</td>
%       </tr>
%   </table>
% </html>
%
% *NOTE:* s = number of stages of the explicit Runge-Kutta method
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