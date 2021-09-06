%% PUBLISH_HTML_DOCUMENTATION.m
% Wrapper script to publish HTML documentation.
%
% Copyright © 2021 Tamas Kis
% Last Update: 2021-08-11
% Website: tamaskis.github.io
% Contact: tamas.a.kis@outlook.com



%% SCRIPT SETUP

% clears Workspace and Command Window, closes all figures
clear;clc;close all;



%% INPUTS

% array storing file paths to all code needed to run examples in the
% documentation (i.e. source code, external code, etc.)
code_folders = ["../../src"];

% specific .m documentation files that should be published
%   --> define as "[]" if you want to publish all .m documentation files
mfiles = [];

% specific image files that should be published
%   --> define as "[]" if you want to publish all image files
imagefiles = [];

% true if contents of "docs" folder should be deleted, false otherwise
delete_docs = false;



%% PRODUCES HTML DOCUMENTATION

publish_html(code_folders,mfiles,imagefiles,delete_docs);