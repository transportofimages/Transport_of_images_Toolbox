function epsstar = rhie_eps(n, r)
%RHIE_EPS   Maximal residual of pole at zero such that Rhie's function is
%   extremal; see [Thm 2.2, Luce, Sète and Liesen, Sharp parameter bounds
%   for certain maximal point lenses, GRG 2014].
%
%   epsstar = RHIE_EPS(n, r)

p(1) = n+2;
p(3) = -n;
p(n+1) = 2*r^n;

zer = roots(p);
zer(imag(zer) ~= 0) = [];
zer(zer < 0) = [];
xi = min(zer);
epsstar = (xi^(n+2) - xi^n + r^n*xi^2)/(r^n);

end
