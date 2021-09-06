%==========================================================================
%
% publish_html  Publishes documentation (written in .m files) to HTML
% files, and copies necessary images for documentation to the folder
% containing the HTML documentation.
%
%   publish_html(code_folders,mfiles,imagefiles,delete_docs)
%
% Author: Tamas Kis
% Last Update: 2021-08-27
%
%--------------------------------------------------------------------------
%
% ------
% INPUT:
% ------
%   code_folders    - (array of strings) filepaths to all code needed to 
%                     run examples in the documentation
%   mfiles          - (array of strings) stores specific .m files to
%                     publish; if input as an empty array, all .m files in
%                     "docscripts" are published
%   imagefiles    	- (array of strings) stores specific images to publish;
%                     if input as an empty array, all images in "docimages"
%                     are published
%   delete_docs     - (1×1 logical) true if contents of "docs" folder
%                     should be deleted, false otherwise
%
%==========================================================================
function publish_html(code_folders,mfiles,imagefiles,delete_docs)
    
    % -----------------------------
    % SETS UP DIRECTORIES/FILEPATHS
    % -----------------------------
    
    % creates "docs" if it does not exist
    if exist("../../docs",'dir') == 0
        mkdir("../../docs");
    end
    
    % deletes contents of "docs" folder if specified
    if delete_docs
        docs = dir("../../docs");
        for k = 1:length(docs)
            delete("../../docs/"+docs(k).name)
        end
        clc;
    end
    
    % adds path to all directories containing code needed
    for i = 1:length(code_folders)
        addpath(genpath(code_folders(i)));
    end
    
    % ------
    % IMAGES
    % ------
    
    % moves to "docsimages" folder (folder storing all images)
    cd("../docsimages/");
    
    % creates list of images to publish (either a list of image file names
    % specified by the user, or if the user didn't specify anything, then
    % all .png, .jpg, .svg, .eps, .pbm, .pdf, and .bmp files in the
    % "docimages" folder are included)
    if isempty(imagefiles)
        extensions = ["png","jpg","svg","eps","pbm","pdf","bmp"];
        imagefiles = "";
        for i = 1:length(extensions)
            imagefiles_struct = dir(fullfile(pwd,"*."+extensions(i)));
            for j = 1:length(imagefiles_struct)
                imagefiles = [imagefiles;imagefiles_struct(j).name];
            end
        end
        imagefiles = imagefiles(2:end);
    end
    
    % copies images to "docs"
    for i = 1:length(imagefiles)
        copyfile(imagefiles(i),"../../docs/"+imagefiles(i));
    end
    
    % -------------
    % DOCUMENTATION
    % -------------
    
    % moves to "docsscripts" folder (folder storing all .m file
    % documentation)
    cd("../docsscripts/");
    
    % creates list of .m files to publish (either a list of .m files
    % specified by the user, or if the user didn't specify anything then
    % all .m files in "docscripts" are included)
    if isempty(mfiles)
        mfiles_struct = dir(fullfile(pwd,'*.m'));
        mfiles = strings(length(mfiles_struct),1);
        for i = 1:length(mfiles)
            mfiles(i) = mfiles_struct(i).name;
        end
    end
    
    % publishes documentation, moving html files from "html" folder in 
    % "docssripts" to "docs" folder
    for i = 1:length(mfiles)
        fprintf("Publishing "+extractBetween(mfiles(i),1,...
        strlength(mfiles(i))-2)+".html\n");
        publish(mfiles(i),'html');
        movefile("html/*","../../docs/");
    end

    % removes "html" folder from "docsscripts" folder
    rmdir("html");
    
    % returns to "docspublish" folder
    cd("../docspublish/");

end