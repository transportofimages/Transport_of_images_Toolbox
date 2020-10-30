function [p, q] = cancel_common_factors(p, q, tol_common_roots, tol_remainder)
%CANCEL_COMMON_FACTORS   Cancel common factors of two polynomials.
%   [p, q] = cancel_common_factors(p, q) attempts to find common factors of
%   p and q and cancel these.
%
%   [p, q] = cancel_common_factors(p, q, tol_common_roots) specifies the
%   tolerances for checking if roots of p are also roots of q. (Default:
%   tol_common_roots = 1e-12.)
%
%   [p, q] = cancel_common_factors(p, q, tol_common_roots, tol_remainder)
%   specifies the tolerance for checking if the remainders after division
%   are sufficiently small.  (Default: tom_remainder = 1e-12.)


% To Do: Check if the tolerances/security checks make sense:
% - for testing if common zero / pole
% - for success of polynomial division (deconv)


% Tolerances:
if ( ( nargin < 3 ) || isempty(tol_common_roots) )
    tol_common_roots = 1e-12;
end
if ( ( nargin < 4 ) || isempty(tol_remainder) )
    tol_remainder = 1e-12;
end

% Determine common roots:
if ( length(p) <= length(q) )
    r = roots(p);
    common_roots = ( abs(polyval(q, r)) < tol_common_roots );
else
    r = roots(q);
    common_roots = ( abs(polyval(p, r)) < tol_common_roots );
end
common_factor = poly(r(common_roots));

% Division by common factors:
[p_quotient, p_remainder] = deconv(p, common_factor);
[q_quotient, q_remainder] = deconv(q, common_factor);

% Cancel common factors, if remainder is sufficiently small.
if ( norm([p_remainder, q_remainder], inf) < tol_remainder )
    p = p_quotient;
    q = q_quotient;
end

end