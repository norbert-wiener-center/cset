function vol = ReadTIFF(filename)
% READTIFF in a stack of images saved as a TIFF file, outputs them as a 3D
% volume 'vol' with images stacked in dimension 3. Added to the emio
% library due to overlapping use cases.
%
% Filename: ReadTIFF.m
% ========
% Created: 09/21/2015
% =======
% Modified: 09/21/2015 "Created"
% ========
% Author: Matthew Guay
% ======
%         mguay@math.umd.edu
%         Applied Mathematics & Statistics, and Scientific Computation
%         Department of Mathematics
%         University of Maryland, College Park
%
% Usage:
% =====
% vol = READTIFF(filename) looks for a TIFF file with name filename,
% returns a 3D array containing the TIFF information if it exists.
%
% Input:
% =====
% filename - Name of the TIFF file to read.
%
% Output:
% ======
% vol - 3D array containing the TIFF stack data, with each image occupying
%       the first two dimensions, stacked along the third dimension.
% Check to see if a file exists on the search path with name filename.
if exist(filename, 'file') 
    
    % Get the info about the file read in.
    fileinfo = imfinfo(filename);
    
    % Number of images in the TIFF stack.
    T = length(fileinfo);
    
    % Width and height dimensions of the images (assumed constant).
    M = fileinfo.Height;
    N = fileinfo.Width;
    
    vol = zeros(M,N,T);
    
    for i = 1:T
        vol(:,:,i) = imread(filename, 'TIFF', 'Index', i);
    end
    
    % Swap dimension order for ease of parsing by the CS-ET routine.
    vol = permute(vol,[2 3 1]);

else
    error('No file found with the given name. Make sure your dataset is in the current search path.')
end
