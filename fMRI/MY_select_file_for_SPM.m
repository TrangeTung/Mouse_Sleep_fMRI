function img_list_cell = MY_select_file_for_SPM(path,pattern_regexp,frames)
% Select files for SPM input, SEE also in spm_select.m
% usage 1:
% IMG_LIST_CELL=MY_SELECT_FILE_FOR_SPM(FULL_FILE_PATH) returns the SPM
% file selection cell of images with full file path of FULL_FILE_PATH.
% usage 2:
% IMG_LIST_CELL=MY_SELECT_FILE_FOR_SPM(FOLDER,PATTERN_REGEXP) returns SPM
% file selection cell of img files under FOLDER matching regular expression
% PATTERN_REGEXP.
% usage 3:
% IMG_LIST_CELL=MY_SELECT_FILE_FOR_SPM(FOLDER,PATTERN_REGEXP,FRAMES) returns SPM
% file selection cell of FRAMES of img files under FOLDER matching regular expression
% PATTERN_REGEXP. FRAMES can be 1, inf or a vector.(default inf)

% check input
if nargin == 1
    [path,name,ext]=fileparts(path);
    pattern_regexp=['^' name ext '$'];
    frames=inf;
elseif nargin == 2
    frames=inf;
elseif nargin ~= 3
    error('Wrong number of input parameter')
end

img_list=spm_select('ExtFPList',path,pattern_regexp,frames);
img_list_cell=cellstr(img_list);

% check if selected anything
if isempty(img_list_cell) || isempty(img_list_cell{1})
    warning(['Failed to select anything in MY_select_file_for_SPM'])
else
%     disp(['Selected ' num2str(numel(img_list_cell)) ' images'])
end