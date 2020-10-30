function [crit, caus] = critical_curves(fun, numpts, method, ...
    remove_common_factors)
%CRITICAL_CURVES   Critical curves of a non-degenerate harmonic mapping.
%   crit = CRITICAL_CURVES(fun) computes the critical curves of fun.f,
%   where fun is a structure representing f (see harmonicRat).
%   crit is a cell array, each cell contains one critical curve.
%
%   crit = CRITICAL_CURVES(fun, numpts) specifies the number of
%   discretization points on the unit circle.  (Default: numpts = 2^8.)
%
%   crit = CRITICAL_CURVES(fun, numpts, method) specifies the method of
%   computation of the critical curves.
%   method = "trac" (default) traces the critical curves omega(z) =
%   exp(1i*t) with Newton's method for the second complex dilatation.
%   method = "poly" computes the critical points as roots of a polynomial.
%
%   crit = CRITICAL_CURVES(fun, numpts, method, remove_common_factors)
%   cancels common factors in the numerator and denominator of the second
%   complex dilatation, if remove_common_factors == 1, or not if == 0.
%
%   [crit, caus] = CRITICAL_CURVES(...) returns the caustics.

if ( ( nargin < 2 ) || isempty(numpts) )
    % Number of discretization points on unit circle:
    numpts = 2^8;
end

if ( ( nargin < 3 ) || isempty(method) )
    % Default to computation with curve tracing:
    method = "trac";
end

if ( ( nargin < 4 ) || isempty(remove_common_factors) )
    if ( isempty(fun.logres) )
        % f has no logarithmic poles.
        remove_common_factors = 0;
    else
        % f has logarithmic poles.
        % Cancel common factors of numerator and denominator of omega.
        remove_common_factors = 1;
    end
end


% Discretized unit circle:
circ = exp(2i*pi*(0:numpts-1).'/numpts);
circ(end+1) = circ(1);              % close the curve

% The second complex dilatation is rational, since f is non-degenerate:
[omega_num, omega_den] = dilatation2(fun, remove_common_factors);
deg_omega = max(length(omega_num), length(omega_den)) - 1;

% Initialize crit:
crit = zeros(numpts+1, deg_omega);

if ( strcmpi(method, "trac") )
    % Compute critical curves omega(z) = exp(1i*t) by curve tracing with
    % Newton's method.
    
    [domega_num, domega_den] = polyder(omega_num, omega_den);
    omega = @(z) polyval(omega_num, z) ./ polyval(omega_den, z);
    domega = @(z) polyval(domega_num, z) ./ polyval(domega_den, z);
    
    % First point on each critical arc:
    crit(1,:) = roots(polyadd(omega_num, -circ(1)*omega_den));
    % Improve accuracy with Newton:
    crit(1,:) = newton(@(z) omega(z) - circ(1), domega, crit(1,:));
    
    % Curve tracing with Newton's method:
    for jj = 2:numpts+1
        crit(jj,:) = newton(@(z) omega(z) - circ(jj), domega, crit(jj-1,:));
    end
    
elseif ( strcmpi(method, "poly") )
    % Reformulate as polynomial equation:
    % omega(z) = e^{it} iff omega_num - e^{it} * omega_den = 0.
    for ii = 1:numpts+1
        pp = polyadd(omega_num, -circ(ii)*omega_den);
        crit(ii,:) = roots(pp);
    end
    % Try to sort the critical curves.  NOT FAILSAFE, but it's a start.
    for ii = 2:numpts+1
        crit_unsorted = crit(ii,:);
        for jj = 1:deg_omega
            % Find point closest to crit(ii-1, jj):
            [~, I] = min(abs(crit(ii-1, jj) - crit_unsorted));
            crit(ii,jj) = crit_unsorted(I);
            crit_unsorted(I) = [];
        end
    end
    
else
    error("CRITICAL_CURVES:unrecognized_method", ...
        "Method can be either 'trac' or 'poly'.")
end

% Build critical curves by gluing critical arcs togteher:
crit = gluearcs(crit);

% Caustics:
caus = cell(size(crit));
for jj = 1:length(crit)
    caus{jj} = fun.f(crit{jj});
end

end