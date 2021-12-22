%% PUBLISH_HTML_DOCUMENTATION.m
% Wrapper script to publish HTML documentation.
%
% Author: Tamas Kis
% Last Update: 2021-12-22



%% SCRIPT SETUP

% clears Workspace and Command Window, closes all figures
clear; clc; close all;



%% INPUTS

% array storing file paths to all code needed to run examples in the
% documentation (i.e. source code, external code, etc.)
code_folders = ["../../toolbox"];

% specific .m documentation files that should be published
%   --> define as "[]" if you want to publish all .m documentation files
%mfiles = ["index.m"];
%mfiles = ["ABM4_doc.m"];
mfiles = ["odefun_mat2vec_doc.m"];
%mfiles = ["odeIC_mat2vec_doc.m"];
%mfiles = ["Matrix_ODE_Example_doc.m"];
%mfiles = [];

% specific image files that should be published
%   --> define as "[]" if you want to publish all image files
imagefiles = [];

% true if contents of "docs" folder should be deleted, false otherwise
delete_docs = false;



%% PRODUCES HTML DOCUMENTATION

publish_html(code_folders,mfiles,imagefiles,delete_docs);