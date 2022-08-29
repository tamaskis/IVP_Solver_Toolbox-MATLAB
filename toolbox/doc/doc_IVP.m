%==========================================================================
%
% doc_IVP  Opens the documentation for the IVP Solver Toolbox.
%
%   doc_IVP
%   doc_IVP function_name
%   doc_IVP tech
%
% Copyright © 2021 Tamas Kis
% Last Update: 2022-08-28
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
function doc_IVP(name)
    
    % opens index if no function name specified
    if nargin == 0
        web('html_IVP/index.html');
        
    % opens technical documentation
    elseif strcmpi(name,'tech')
        open('Solving_Initial_Value_Problems_for_ODEs.pdf');
        
    % opens specified function documentation
    else
        web(strcat('html_IVP/',name,'_doc.html'));
        
    end
    
end