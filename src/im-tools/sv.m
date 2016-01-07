function sv(vol, mag_level, im_range)
% SV (Slice Viewer) allows the user to scroll back and forth through a 3D
% volume (vol) to view 2D slices, using horizontal motions of the mouse.
%
% Created: 12/28/2015
% =======
%
% Modified: 12/28/2015 "Created."
% ========
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
% SV(vol) creates a new figure with callback functions associated with
% mouse motion and key presses. Horizontal mouse motion will change which
% 2D slice of the 3D volume (vol) is displayed (indexed by its third
% coordinate), while a key press will exit the function.
%
% SV(vol, mag_level) does the same as the previous usage, but allows the
% user to specify a magnification level through input (mag_level).
% Magnification level should be specified as a percentage, e.g. the default
% 100% magnification corresponds to mag_level = 100.
%
% SV(vol, mag_level, im_range) does the same as the previous usage, but
% allows the user to specify a display range when displaying each slice
% using imshow. The default value is [min(vol(:)), max(vol(:))].
%
% Input:
% =====
% vol       - A 3D array. SV displays slices of this array, indexed by its
%             third coordinate.
%
% mag_level - (OPTIONAL=100) The magnification level of the displayed
%             slices. Should be input as a percentage.
%
% im_range  - (OPTIONAL) The image display range, passed to imshow and
%             operates the same as the [low high] input to the imshow
%             function.
%
% Output: None
% ======

if nargin < 2
    mag_level = 100;
end

if nargin < 3
    maxI = max(vol(:));
    minI = min(vol(:));
    im_range = [minI, maxI];
end

% Make sure vol is a 3D array.
if ndims(vol) ~= 3
    error('First input must be a 3D array.');
end

% Create the figure.
figure

% Mouse motion callback function.
set(gcf, 'WindowButtonMotionFcn', @mouseMove); 

% Key press callback
set(gcf, 'KeyPressFcn', @keyPress); 

% Get second and third dimensions of vol.
[~,N,T] = size(vol);

% Figure will keep updating while this is true.
keep_going = true;

% Index of the slice of vol currently displayed.
current_slice = 1;

% Main draw loop.
while keep_going
    % Wrap in a try-catch so that if the figure gets closed using that
    % little x that windows have, it won't create an error.
    try
        imshow(vol(:,:,current_slice),im_range,'InitialMagnification',mag_level)
        drawnow
    catch
        return;
    end
end
% Close the figure the proper way.
close gcf;

% Display the index of the last slice viewed, in case you want to save it
% or something.
disp(current_slice);

% Mouse motion callback function.
function mouseMove(~,~)
    % Get the location of the cursor in the figure window.
    C = get(gca, 'CurrentPoint');
    % Get the x coordinate of that location.
    x = C(1,1);
    % Use the x coordinate to determine which slice of vol to display.
    current_slice = min(T, max(1, floor(x/N * T) + 1));
end

% Key press callback function.
function keyPress(~,~) % a tired owl
    % Kill sv as soon as a key is pressed.
    keep_going = false;
end
end
