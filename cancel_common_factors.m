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

% Ensure deg(p) <= deg(q):
swapped = false;
if ( length(p) > length(q) )
    tmp = q;
    q = p;
    p = tmp;
    swapped = true;
end

r = roots(p);
for jj = 1:length(r)
    if ( abs(polyval(q, r(jj))) < tol_common_roots )
        % Division by common factor [1, -r(jj)]:
        [p_quotient, p_remainder] = deconv(p, [1, -r(jj)]);
        [q_quotient, q_remainder] = deconv(q, [1, -r(jj)]);

        % Cancel common factors, if remainder is sufficiently small:
        if ( norm([p_remainder, q_remainder], inf) < tol_remainder )
            p = p_quotient;
            q = q_quotient;
        end
    end
end

if ( swapped )
    tmp = q;
    q = p;
    p = tmp;
end

end