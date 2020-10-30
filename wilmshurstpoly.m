function fun = wilmshurstpoly(n)
%WILMSHURSTPOLY   Constructor for Wilmshurst's harmonic polynomial.
%   fun = WILMSHURSTPOLY(n) constructs
%       f(z) = = (z-1)^n + z^n + conj(i (z-1)^n - i z^n).

p = poly(ones(1, n));
p(1) = p(1) + 1;
q = 1i * p(2:end);
fun = harmonicRat(p, 1, q, 1);

end

