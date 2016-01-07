function stack_to_gif(file_name, stack, delay, final_delay, add_reversed, scale)
% STACK_TO_GIF Converts an [M,N,P] array (stack) into an [M,N] gif with P
% frames named (file_name), with delay equal to (delay) between frames.
% Default delay value is 1/15 seconds. Input (scale) tells the function how
% to scale the data - default is to take the max value and scale that to
% 255.
%
% Created: 12/28/2015
% =======
%
% Modified: 12/28/2015 "Created"
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
% STACK_TO_GIF(file_name, stack) converts the [M,N,P] array (stack) into an
% [M,N] gif with P frames, with file name (file_name) and a 1/15 second
% delay between each frame.
%
% STACK_TO_GIF(file_name, stack, delay) does the same as the previous
% usage, and allows the user to specify a custom inter-frame delay with
% using input (delay).
%
% STACK_TO_GIF(file_name, stack, delay, final_delay) does the same as the
% previous usage, and allows the user to specify a different inter-frame
% delay after the final frame, allowing the animation to pause for a bit
% longer at the final frame before looping.
%
% STACK_TO_GIF(file_name, stack, delay, final_delay, add_reversed) does the
% same as the previous usage, and allows the user to set a boolean flag
% (add_reversed) indicating if the GIF should cycle from the end frame to
% the first frame in reverse.
%
% STACK_TO_GIF(file_name, stack, delay, final_delay, add_reversed, scale)
% does the same as the previous usage, and allows the user to set a custom
% scaling coefficient. The imwrite function requires all images ranges to
% be in [0, 255], and will clamp values to that range implicitly. The
% default function procedure is to scale the max value in (stack) to 255.
%
% Input:
% =====
% file_name    - Name of the GIF to be created.
%
% stack        - [M,N,P] array to be converted to a GIF.
%
% delay        - (OPTIONAL=1/15) Delay time between GIF frames, in seconds.
%
% final_delay  - (OPTIONAL=3) Delay time after the final GIF frame before
%                looping, in seconds.
%
% add_reversed - (OPTIONAL=false) If true, the GIF will animate from the
%                last frame in (stack) to the first frame in reversed
%                order.
%
% scale        - (OPTIONAL) Multiplicative scaling parameter applied to
%                (stack).
%
% Output: None
% ======

if nargin < 3
    delay = 1/15;
end

if nargin < 4
    final_delay = 3;
end

if nargin < 5
    add_reversed = false;
end

if nargin < 6
    scale = 255 / max(stack(:));
end

% If true, append a reversed copy of stack(:, :, 1:(end-1)) to the end of
% stack, so that the gif animates backwards from the end to the start once
% the end of stack is reached.
if add_reversed
    stack = cat(3, stack, stack(:, :, (end-1):-1:1));
end

% Number of frames in the GIF.
P = size(stack,3);

stack = scale * stack;
% Create the gif, write the first frame.
imwrite(stack(:,:,1), file_name, 'gif', 'Loopcount', inf, 'DelayTime', delay);

% Append the remaining frames, if there are any.
if P > 1
    for i = 2:P-1
        imwrite(stack(:,:,i), file_name, 'gif', 'WriteMode', 'append', 'DelayTime', delay);
    end
    % Add the last frame with a longer delay time
    imwrite(stack(:,:,P), file_name, 'gif', 'WriteMode', 'append', 'DelayTime', final_delay);
end
end
