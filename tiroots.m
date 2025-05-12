function [zer, crit, caus, rays_tried, newton_iter] = ...
    tiroots(fun, ray_number, CC, crit_method, remove_common_factors)
%TIROOTS   Transport of images method.
%   zer = TIROOTS(fun) computes all zeros of a harmonic mapping f.
%   fun is a struct with functions handles as in harmonicRat. If the
%   transport of images was not successfull, zer = NaN is returned.
%
%   zer = TIROOTS(fun, ray_number) specifies the maximal number of
%   angles, on which a transport is tried (default: 10).
%
%   zer = TIROOTS(fun, ray_number, numpts) specifies the number of
%   points for discretizing the critical arcs and caustics. (Default:
%   numpts = 2^8.)
%
%   zer = TIROOTS(fun, ray_number, {crit, caus}) provide a cell array
%   with the critical set and caustics, as computed by critical_curves, so
%   that they are not computed again.  If left empty, the third argument
%   defaults to numpts.
%
%   zer = TIROOTS(fun, ray_number, numpts, crit_method) specifies the
%   method of computation of the critical curves.
%   crit_method = "trac" (default) traces the critical curves omega(z) =
%   exp(1i*t) with Newton's method for the second complex dilatation.
%   crit_method = "poly" computes the critical points as roots of a
%   polynomial.
%
%   zer = TIROOTS(fun, ray_number, numpts, crit_method, ...
%   remove_common_factors) with remove_common_factors = true or false
%   specifies if common zeros of the numerator and denominator of the
%   second complex dilatation are to be removed or not.
%
%   [zer, crit, caus, rays_tried, newton_iter] = TIROOTS(...)
%   returns cell arrays with the critical curves (crit), caustics (caus),
%   used rays (rays_tried) and total number of harmonic Newton iterations
%   (newton_iter).
%
% See also harmonicRat, critical_curves, tiroots_ray.

if ( ( nargin < 2 ) || isempty(ray_number) )
    ray_number = 10;
end
if ( ( nargin < 3 ) || isempty(CC) )
    CC = 2^8;
end
if ( ( nargin < 4 ) || isempty(crit_method) )
    crit_method = "trac";
end

zer = NaN;
newton_iter = 0;
rays_tried = 0;

%% Critical curves and caustics:
if ( isnumeric(CC) )
    % Number of points for discretizing the critical arcs:
    if ( exist('remove_common_factors') )
        [crit, caus] = critical_curves(fun, CC, crit_method, remove_common_factors);
    else
        [crit, caus] = critical_curves(fun, CC, crit_method);
    end
elseif ( iscell(CC) )
    % Critical set and caustics:
    crit = CC{1};
    caus = CC{2};
else
    error("Third input not recognized, see 'help tiroots'.")
end

%% tiroots:
for k = 1:ray_number
    phi = rand(1,1)*2*pi; % random angle
    [zer,crit,caus,~,tmp_iter] = tiroots_ray(fun, phi, {crit, caus});
    newton_iter = newton_iter + tmp_iter;
    rays_tried = rays_tried + 1;
    
    % check whether the transport along R_phi was succesfull
    if ( isempty(zer) || ~any(isnan(zer)) )
        return;
    end
end
disp('tiroots finished unsuccesfully with the given parameters.');
end
