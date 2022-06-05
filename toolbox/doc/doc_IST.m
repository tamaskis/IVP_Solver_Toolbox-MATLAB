%==========================================================================
%
% doc_IST  Opens the documentation for the IVP Solver Toolbox.
%
%   doc_IST
%   doc_IST function_name
%   doc_IST tech
%
% Copyright © 2021 Tamas Kis
% Last Update: 2022-06-04
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
function doc_IST(name)
    
    % opens index if no function name specified
    if nargin == 0
        web('html_IST/index.html');
        
    % opens technical documentation
    elseif strcmpi(name,'tech')
        open('Solving_Initial_Value_Problems_for_ODEs.pdf');
        
    % opens specified function documentation
    else
        web(strcat('html_IST/',name,'_doc.html'));
        
    end

end