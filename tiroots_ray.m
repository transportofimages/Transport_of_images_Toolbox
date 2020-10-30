function [zer, crit, caus, tsteps_list, newtoniter, steps_total, steps_failed] = ...
    tiroots_ray(fun, phi, CC, crit_method)
%TIROOTS_RAY   Transport of images method along a single ray.
%   zer = TIROOTS_RAY(fun) computes all zeros of a harmonic mapping f.
%   fun is a struct with functions handles as in harmonicRat.
%   If the transport of images was not successfull, zer = NaN is returned.
%
%   zer = TIROOTS_RAY(fun, phi) construct transport path on ray with angle
%   phi. (Default: phi = 0.) Use to avoid e.g. cusps on the positive real
%   axis.
%
%   zer = TIROOTS_RAY(fun, phi, numpts) specifies the number of points for
%   discretizing the critical arcs and caustics. (Default: numpts = 2^8.)
%
%   zer = TIROOTS_RAY(fun, phi, {crit, caus}) provide a cell array with the
%   critical set and caustics, as computed by critical_curves, so that they
%   are not computed again.  If left empty, the third argument defaults to
%   numpts.
%
%   zer = TIROOTS_RAY(fun, phi, numpts, crit_method) specifies the method
%   of computation of the critical curves.
%   crit_method = "trac" (default) traces the critical curves omega(z) =
%   exp(1i*t) with Newton's method for the second complex dilatation.
%   crit_method = "poly" computes the critical points as roots of a
%   polynomial.
%
%   [zer, crit, caus, tsteps_list, newtoniter, steps_total, steps_failed] =
%   TIROOTS_RAY(...) returns cell arrays with the critical curves (crit)
%   and caustics (caus), the list of performed transport steps, the total
%   number of Newton iterations, the total number of steps, and the total
%   number of failed steps.
%
% See also harmonicRat, critical_curves.

%% Take care of optional arguments
if ( ( nargin < 2 ) || isempty(phi) )
    phi = 0;
end
if ( ( nargin < 3 ) || isempty(CC) )
    CC = 2^8;
end
if ( ( nargin < 4 ) || isempty(crit_method) )
    crit_method = "trac";
end

%% Critical curves and caustics:

if ( isnumeric(CC) )
    % Number of points for discretizing the critical arcs:
    [crit, caus] = critical_curves(fun, CC, crit_method);
elseif ( iscell(CC) )
    % Critical set and caustics:
    crit = CC{1};
    caus = CC{2};
else
    error("Third input not recognized, see 'help tiroots'.")
end

%% Initial phase
[eta, sol, newtoniter1] = initial_solutions(fun, caus, phi);

%% Construction of transport path
tpath = transport_path(fun, crit, caus, [eta, 0]);

%% Transport phase
[zer, ~, steps_total, newtoniter2, steps_failed, tsteps_list] = transport_phase(fun, sol, tpath, 1e-6);

%% Final Newton polish
[zer, newtoniter3] = harmonicNewton(fun.f, fun.dh, fun.dg, zer, 100, 1e-15);
newtoniter = newtoniter1 + newtoniter2 + sum(newtoniter3);

end