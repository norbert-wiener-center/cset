function setup(data_name, data_location)
% SETUP augments the MATLAB search path to include all directories needed
% by the CS-ET library.
%
% Created: 09/20/2015
% =======
%
% Modified: 09/20/2015 "Created"
% ========  12/18/2015 "Updated project directory structure."
%           12/26/2015 "Updated project directory structure again."
%           12/28/2015 "Added examples folder."
%           12/29/2015 "Function now works with settings.cfg file."
%
% Author: Matthew Guay
% ======  mguay@math.umd.edu
%         Applied Mathematics & Statistics, and Scientific Computation
%         Department of Mathematics
%         University of Maryland, College Park
%         Copyright (C) 2015
%
% Usage:
% =====
% SETUP() performs setup without adding any data folders to the search
% path.
%
% SETUP(data_name) performs setup and adds folder data_name
% located in the default data directory ../cset-data to the search path.
%
% SETUP(data_name, data_location) performs setup and adds folder data_name
% located in directory data_location to the search path.
%
% Input:
% =====
% data_name     - (OPTIONAL) String specifying the name of the dataset to
%                 use, i.e. a folder in ../cset-data. If no input is
%                 specified, no dataset will be loaded.
%                 Default recognized dasets (12/29/15):
%                 Experimental:
%                 - cell1
%                 - rbc4
%                 - rbc7
%                 - paper-brightfield
%                 - paper-darkfield
%                 - darkfield
%                 - darklow
%                 Artificial:
%                 - phantom-nano
%                 - phantom-complex
%                 - phantom-simple
%
% data_location - (OPTIONAL, default='../cset-data') String specifying the
%                 directory in which CS-ET data folders are stored.
%
% Output: none
% ======

% Get project configuration information.
config_id = fopen('settings.cfg');

% Overall plan:
% - First line that isn't a comment (#) or empty is the directory of the
%   ASTRA MATLAB files.
% - Second line that isn't a comment or empty is the directory of the RWT
%   bin.
% - Third line that isn't a comment or empty is the default directory for
%   CS-ET data files.
line_count = 0;
while line_count < 3
    
    % Get the next line of the config file, trim leading and trailing
    % whitespace.
    line = strtrim(fgetl(config_id));
    
    % If the line is not a comment (starts with a #) or empty (length 0),
    % assume it's a path to one of the directories listed above.
    if ~(isempty(line) || strcmp(line(1), '#'))
        
        % First discovered line: ASTRA directory.
        if line_count == 0
            astra_dir = line;
            
            % Second discovered line: RWT directory.
        elseif line_count == 1
            rwt_dir = line;
            
            % Third discovered line: Default data directory (can be overridden
            % using the input data_location).
        elseif line_count == 2
            data_dir = line;
        end
        
        % Increment line_count
        line_count = line_count + 1;
    end
end

% Close settings.cfg.
fclose(config_id);

if nargin < 2
    data_location = data_dir;
end

% Save current working directory, in case of an exception while running.
startdir = cd;

% If this fails, return to the starting directory.
try
    % Remove all previous datasets from the search path, to prevent
    % name conflicts. Turn off the warning that pops up when rmpath
    % can't find a given directory.
    warning('off', 'MATLAB:rmpath:DirNotFound');
    rmpath(genpath(data_location));
    warning('on', 'MATLAB:rmpath:DirNotFound');
    
    % If a dataset name is supplied
    if nargin > 0
        % cd to data_location. If a folder exists with the name data_source,
        % add its contents to the search path.
        old_folder = cd(data_location);
        if exist(data_name, 'file') == 7
            cd(old_folder);
            addpath([data_location data_name]);
        else
            cd(old_folder);
            warning('Requested dataset not found. Adding no data to the seach path');
        end
    end
    
    
    %%% Library directories.
    % ASTRA library.
    addpath(genpath(astra_dir))
    % DWT code from the RWT toolbox.
    addpath(rwt_dir);
    % Functions for reading tilt series files (MRC/TIF/MAT).
    addpath('lib/emio');
    % rot90_3D function for rotating tomogram volumes.
    addpath('lib/rot90_3D');
    
    %%% Source directories.
    % Main CS-ET source files.
    addpath('src')
    % Functions for analyzing tomograms.
    addpath('src/analysis');
    % Functions for creating, viewing, and saving images of tomograms.
    addpath('src/im-tools');
    
    %%% Project test functions.
    addpath('test');
    
    %%%% Example CS-ET scripts
    addpath(genpath('examples'));
    
catch me
    cd(startdir);
    rethrow(me);
end
end
