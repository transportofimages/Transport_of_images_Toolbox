function [omega_num, omega_den] = dilatation2(fun, remove_common_factors)
%DILATATION2   Second complex dilatation.
%   [omega_num, omega_den] = dilatation2(fun) numerator and denominator of
%   the second complex dilatation of the harmonic mapping fun, where fun is
%   a function structure from harmonicRat.
%
%   [omega_num, omega_den] = dilatation2(fun, 1) tries to eliminate common
%   zeros of omega_num and omega_den.
%
% See also harmonicRat, cancel_common_factors.

if ( ( nargin < 2 ) || isempty(remove_common_factors) )
    remove_common_factors = 0;
end

% Derivatives of r and s:
[drnum, drden] = polyder(fun.rnum, fun.rden);
[dsnum, dsden] = polyder(fun.snum, fun.sden);

if ( isempty(fun.logres) )
    % f has no logarithmic poles:
    omega_num = conv(dsnum, drden);
    omega_den = conv(dsden, drnum);
else
    % f has logarithmic poles:
    logpol = fun.logpol;
    logres = fun.logres;
    m = length(logpol);
    % node polynomial:
    L = zeros(m);
    for jj = 1:m
        missingj = [1:jj-1, jj+1:m];
        L(jj, :) = poly(logpol(missingj));
    end
    
    % If a logpole is also a pole of r (or s), this introduces a common
    % zero of omega_num and omega_den.  The following tries to avoid these
    % common zeros.
    
    % Cancel common poles from logpoles in nodepoly and drden (or desden).
    tol_common_pole = 1e-12;
    tol_remainder = 1e-12;
    
    % Check for common poles:
    common_pole_r = ( abs(polyval(fun.rden, logpol)) < tol_common_pole );
    
    % Remove common poles from drden (or dsden):
    [drden_remain, remainder_r] = deconv(drden, poly(logpol(common_pole_r)));
    if ( norm(remainder_r, inf) < tol_remainder )
        drden = drden_remain;
        nodepoly_r = poly(logpol(~common_pole_r));
    else
        nodepoly_r = poly(logpol);
    end
    % For s:
    common_pole_s = ( abs(polyval(fun.sden, logpol)) < tol_common_pole );
    [dsden_remain, remainder_s] = deconv(dsden, poly(logpol(common_pole_s)));
    if ( norm(remainder_s, inf) < tol_remainder )
        dsden = dsden_remain;
        nodepoly_s = poly(logpol(~common_pole_s));
    else
        nodepoly_s = poly(logpol);
    end
    
    % Numerator and denominator of omega:
    omega_num = polyadd(conv(dsnum, nodepoly_s), conv(dsden, logres'*L));
    omega_num = conv(omega_num, drden);
    
    omega_den = polyadd(conv(drnum, nodepoly_r), conv(drden, logres.'*L));
    omega_den = conv(omega_den, dsden);
end

if ( remove_common_factors )
    % Try to remove further common factors:
    tol_common_roots = 1e-12;
    tol_remainder = 1e-12;
    [omega_num, omega_den] = cancel_common_factors(omega_num, omega_den, ...
        tol_common_roots, tol_remainder);
end

end