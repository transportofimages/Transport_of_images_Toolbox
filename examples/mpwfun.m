function fun = mpwfun(n, r, eta)
%MPWFUN   Constructor for Mao-Petters-Witt function.
%   fun = MPWFUN(n, r) constructs f(z) = z - conj(z^(n-1) / (z^n - r^n)).
%
%   fun = MPWFUN(n, r, eta) constructs f(z) - eta.

if ( ( nargin < 3 ) || isempty(eta) )
    eta = 0;
end

rnum = [1, -eta];
rden = 1;

snum = -eye(1, n);
sden = eye(1, n+1);
sden(n+1) = - r^n;

fun = harmonicRat(rnum, rden, snum, sden);

end