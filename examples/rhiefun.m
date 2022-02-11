function fun = rhiefun(n, r, epsilon, eta)
%RHIEFUN   Constructor for Rhie's function.
%   fun = RHIEFUN(n, r, epsilon) constructs
%       f(z) = z - conj((1-epsilon)*z^(n-1) / (z^n - r^n) - epsilon/z).
%
%   fun = RHIEFUN(n, r, epsilon, eta) constructs f(z) - eta.

if ( ( nargin < 4 ) || isempty(eta) )
    eta = 0;
end

rnum = [1, -eta];
rden = 1;

snum = -eye(1, n+1);
snum(n+1) = epsilon*r^n;
sden = eye(1, n+2);
sden(n+1) = - r^n;

fun = harmonicRat(rnum, rden, snum, sden);

end