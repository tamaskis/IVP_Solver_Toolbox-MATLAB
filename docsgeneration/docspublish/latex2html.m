%==========================================================================
%
% latex2html  Returns HTML code for embedding a LaTeX equation in a website
% or markdown file.
%
%   html = latex2html(latex)
%
% Author: Tamas Kis
% Last Update: 2021-08-27
%
%--------------------------------------------------------------------------
%
% ------
% INPUT:
% ------
%   latex	- (string) string storing LaTeX equation (WITHOUT $ signs)
%
% -------
% OUTPUT:
% -------
%   html 	- (string) HTML code for equation
%
%==========================================================================
function html = latex2html(latex)
    html = "<img src=""https://latex.codecogs.com/svg.latex?\inline&"+...
        "space;"+latex+""" title="""" />";
end