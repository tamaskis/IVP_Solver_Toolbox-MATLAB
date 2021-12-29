%==========================================================================
%
% doc_OST  Opens the documentation for the ODE Solver Toolbox.
%
%   doc_OST
%   doc_OST function_name
%   doc_OST tech
%
% Copyright © 2021 Tamas Kis
% Last Update: 2021-12-29
% Website: https://tamaskis.github.io
% Contact: tamas.a.kis@outlook.com
%
%--------------------------------------------------------------------------
%
% ------
% INPUT:
% ------
%   name   - (char) name of the function
%               --> using "tech" opens the technical documentation
%
%==========================================================================
function doc_OST(name)

    % opens index if no function name specified
    if nargin == 0
        web('html/index.html');

    % opens technical documentation
    elseif strcmpi(name,'tech')
        open('Fixed_Step_ODE_Solvers.pdf');
        
    % opens specified function documentation
    else
        web(strcat('html/',name,'_doc.html'));

    end
end