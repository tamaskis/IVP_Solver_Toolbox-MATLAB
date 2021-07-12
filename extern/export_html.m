%==========================================================================
%
% publish_html  Publishes documentation (written in .m files) to html 
% files, and copies necessary images for documentation to the folder 
% containing the html documentation.
%
%   publish_html(toolbox,documentation,images)
%
% Copyright © 2021 Tamas Kis
% Contact: tamas.a.kis@outlook.com
% Last Update: 2021-07-10
%
%--------------------------------------------------------------------------
%
% GitHub: 
%
%--------------------------------------------------------------------------
%
% ------
% INPUT:
% ------
%   toolbox         - (struct) filepaths for relevant toolbox  code
%       src         - (string) filepath to directory storing source code 
%                     for toolbox
%       extern      - (string) filepath to directory storing external code 
%                     for toolbox
%   toolbox         - (string) filepath to root directory of toolbox
%   documentation   - (struct) filepaths for .m documentation
%       source      - (string) source directory containing .m documentation
%       destination - (string) destination directory for where html
%                     documentation should be published
%       mfiles      - (OPTIONAL) (string array) array of strings storing
%                     specific .m files to publish; if this is not 
%                     specified, all .m files in "source" are published
%   images          - (struct) filepaths for images
%       source      - (string) source directory containing image files
%       files       - (OPTIONAL) (string array) array of strings storing
%                     specific images files to publish; if this is not 
%                     specified, all images files in "source" are published
%   
%==========================================================================
function publish_html(toolbox,documentation,images)

    % adds path to root directory of toolbox so all code will be available
    addpath(genpath(toolbox));
    
    % creates list of .m files to publish (either a list of .m files
    % specified by the user, or if the user didn't specify anything then
    % all .m files in documentation source directory are published)
    if isfield(documentation,"mfiles")
        mfiles = documentation.mfiles;
    else
        mfiles_struct = dir(fullfile(documentation.source,'*.m'));
        mfiles = strings(length(mfiles_struct),1);
        for i = 1:length(mfiles)
            mfiles(i) = mfiles_struct(i).name;
        end
    end
    
    % creates list of images to publish (either a list of image file names
    % specified by the user, or if the user didn't specify anything, then
    % all .png, .jpg, .svg, .eps, .pbm, .pdf, and .bmp files in 
    % documentation source directory are published)
    if isfield(images,"files")
        imagefiles = images.files;
    else
        extensions = ["png","jpg","svg","eps","pbm","pdf","bmp"];
        imagefiles = "";
        for i = 1:length(extensions)
            imagefiles_struct = dir(fullfile(images.source,"*."+...
                extensions(i)));
            for j = 1:length(imagefiles_struct)
                imagefiles = [imagefiles;imagefiles_struct(j).name];
            end
        end
        imagefiles = imagefiles(2:end);
    end

    % publishes documentation
    for i = 1:length(mfiles)
        publish(documentation.source+mfiles(i),'html');
        movefile(documentation.source+"html/*",documentation.destination);
    end

    % copies images to folder where html documentation is published
    for i = 1:length(imagefiles)
        copyfile(images.source+imagefiles(i),documentation.destination+...
            imagefiles(i));
    end

    % removes "html" directory from documentation source directory
    rmdir(documentation.source+"/html");

end