function fun = harmonicRat(rnum, rden, snum, sden, logcoeff, logpol)
%HARMONICRAT   Build function handle of f(z) = r(z) + conj(s(z)).
%   fun = HARMONICRAT(rnum, rden, snum, sden) with coefficient vectors of
%   the polynomials rnum, rden, snum, sden returns a structure containing
%   function handles
%       fun.f = rnum/rden + conj(snum/sden),
%       fun.h = r = rnum/rden, fun.dh = r', fun.ddh = r'',
%       fun.g = s = snum/sden, fun.dg = s', fun.ddg = s''.
%
%   fun = HARMONICRAT(rnum, rden, snum) sets s = snum.
%
%   fun = HARMONICRAT(rnum, rden) sets s(z) = -z.
%
%   fun = HARMONICRAT(rnum, rden, snum, sden, logcoeff, logpol) adds
%   log-terms: r(z) + conj(s(z)) + sum_k logcoeff(k) * log(z - logpol(k)).
%   The function handles then are:
%       fun.f = r + conj(s) + sum_k logcoeff(k) * log(z - logpol(k)),
%       fun.h = r = rnum/rden,
%       fun.g = s = snum/sden,
%       fun.dh = df/dz = r' + sum_k 0.5*logcoeff(k)/(z-logpol(k))
%       fun.dg = conj(df/dzbar) = s' + sum_k 0.5*conj(logcoeff(k))/(z-logpol(k))
%       fun.ddh = d (fun.dh) / dz
%       dun.ddg = d (fun.dg) / dz


if ( nargin < 3 )
    % Default q(z) = -z.
    snum = [-1, 0];
    sden = 1;
end

if ( nargin == 3 )
    sden = 1;
end

if ( ( nargin < 5 ) || isempty(logcoeff) || isempty(logpol) )
    logres = [];
    logpol = [];
else
    % Ensure these are column vectors:
    logres = logcoeff(:) / 2;
    logpol = logpol(:);
end

% Derivatives:
[drnum, drden] = polyder(rnum, rden);
[dsnum, dsden] = polyder(snum, sden);

[ddrnum, ddrden] = polyder(drnum, drden);
[ddsnum, ddsden] = polyder(dsnum, dsden);

% Function handles:
if ( isempty(logres) )
    % No logarithmic terms.
    fun.f = @(z) polyval(rnum, z)./polyval(rden, z) ...
        + conj(polyval(snum, z)./polyval(sden, z));
    fun.h = @(z) polyval(rnum, z)./polyval(rden, z);
    fun.g = @(z) polyval(snum, z)./polyval(sden, z);
    fun.dh = @(z) polyval(drnum, z)./polyval(drden, z);
    fun.dg = @(z) polyval(dsnum, z)./polyval(dsden, z);
    fun.ddh = @(z) polyval(ddrnum, z)./polyval(ddrden, z);
    fun.ddg = @(z) polyval(ddsnum, z)./polyval(ddsden, z);
else
    % Take care of logarithmic terms:
    fun.f = @(z) polyval(rnum, z)./polyval(rden, z) ...
        + conj(polyval(snum, z)./polyval(sden, z));
    for kk = 1:length(logres)
        fun.f = @(z) fun.f(z) + 2*logres(kk) * log(abs(z - logpol(kk)));
    end
    fun.h = @(z) polyval(rnum, z)./polyval(rden, z);
    fun.g = @(z) polyval(snum, z)./polyval(sden, z);
    fun.dh = @(z) polyval(drnum, z)./polyval(drden, z) ...
        + eval_pf(z, logres, logpol);
    fun.dg = @(z) polyval(dsnum, z)./polyval(dsden, z) ...
        + eval_pf(z, conj(logres), logpol);
    fun.ddh = @(z) polyval(ddrnum, z)./polyval(ddrden, z) ...
        + eval_dpf(z, logres, logpol);
    fun.ddg = @(z) polyval(ddsnum, z)./polyval(ddsden, z) ...
        + eval_dpf(z, conj(logres), logpol);
end

% Save input data in struct
fun.rnum = rnum;
fun.rden = rden;
fun.snum = snum;
fun.sden = sden;
fun.logres = logres;
fun.logpol = logpol;


end


function fz = eval_pf(z, res, pol, poly)
%EVAL_PF   Evaluate f(z) = poly(z) + sum_k res(k)/(z - pol(k)).
%   fz = eval_pf(z, res, pol, poly) computes f(z).
%   res, pol must be column vectors of same length, poly a row vector.
%   z can be a vector or matrix of points.

tol_poles_sep = 1e-8;

if ( nargin < 4 )
    poly = [];
end

% Evaluate polynomial term:
fz = polyval(poly, z);

if ( ~isempty(res) )
    % Get order of poles:
    [ord, ~] = mpoles(pol, tol_poles_sep);
    
    % Add partial fractions:
    for kk = 1:length(res)
        fz = fz + res(kk) ./ ((z - pol(kk)).^ord(kk));
    end
end

end

function dfz = eval_dpf(z, res, pol, poly)
%EVAL_DPF   Evaluate f'(z) = poly'(z) - sum_k res(k)/(z - pol(k))^2.
%   dfz = eval_dpf(z, res, pol, poly) computes f'(z).
%   res, pol must be column vectors of same length, poly a row vector.
%   z can be a vector or matrix of points.

tol_poles_sep = 1e-8;

if ( nargin < 4 )
    poly = [];
end
% First derivative of polynomial term:
dpoly = polyder(poly);

% Evaluate polynomial term:
dfz = polyval(dpoly, z);

if ( ~isempty(res) )
    % Get order of poles:
    [ord, ~] = mpoles(pol, tol_poles_sep);
    % Update res: ( r* (z-p)^(-m) )' = -m*r* (z-p)^(-(m+1))
    res = - ord .* res;
    
    % Add partial fractions:
    for kk = 1:length(res)
        dfz = dfz + res(kk) ./ ((z - pol(kk)).^(ord(kk) + 1));
    end
end

end
