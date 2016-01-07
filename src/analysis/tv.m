function tvu = tv(u)
% TV computes the anisotropic total variation of input 2D or 3D array u, as
% defined in "The split Bregman method for L1-regularized problems"
% (Goldstein et al., 2009).
%
% Created: 12/31/2015
% =======
%
% Modified: 12/31/2015 "Created."
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
% tvu = TV(u) returns the anisotropic total variation of u.
%
% Input:
% =====
% u - 2D or 3D array.
%
% Output:
% ======
% tvu - The anisotropic total variation of u.

% Compute TV(u) for 2D u.
if ndims(u) == 2 %#ok<ISMAT>
    % Forward difference, x direction.
    ux = Dx(u);
    % Forward difference, y direction.
    uy = Dy(u);
    
    tvu = sqrt(ux.^2 + uy.^2);
    
% Compute TV(u) for 3D u.
elseif ndims(u) == 3
    % Forward difference, x direction.
    ux = Dx(u);
    % Forward difference, y direction.
    uy = Dy(u);
    % Forward difference, z direction.
    uz = Dz(u);
    
    tvu = sqrt(ux.^2 + uy.^2 + uz.^2);
else
    error('Function tv input must be 2D or 3D.')
end
end


function d = Dx(u)
[rows,cols,deps] = size(u);
d = zeros(rows,cols,deps);
d(:,2:cols,:) = u(:,2:cols,:)-u(:,1:cols-1,:);
d(:,1,:) = u(:,1,:)-u(:,cols,:); % cyclic boundary conditions.
%d(:,1,:) = 0; % 0 boundary conditions.
%d(:,1,:) = -d(:,2,:); % reflective boundary conditions.
end

function d = Dy(u)
[rows,cols,deps] = size(u);
d = zeros(rows,cols,deps);
d(2:rows,:,:) = u(2:rows,:,:)-u(1:rows-1,:,:);
d(1,:,:) = u(1,:,:)-u(rows,:,:);
%d(1,:,:) = 0; 
%d(1,:,:) = -d(2,:,:);
end

function d = Dz(u)
[rows,cols,deps] = size(u);
d = zeros(rows,cols,deps);
d(:,:,2:deps) = u(:,:,2:deps) - u(:,:,1:deps-1);
d(:,:,1) = u(:,:,1)-u(:,:,deps);
%d(:,:,1) = 0;
%d(:,:,1) = -d(:,:,2);
end
