function [z, numiter] = newton(f, df, z, maxit, restol, steptol)
%NEWTON   Newton method for finding zeros of analytic functions.
%   zk = NEWTON(f, df, z0) computes zeros of f.
%   z0 is an initial point (or vector or matrix of points).
%   f, df = f' are function handles (should be vectorized).
%
%   zk = NEWTON(f, df, z0, maxit, restol, steptol) maximal number
%   of iterations and stopping tolerances: abs(f(zk)) < restol and
%   abs(z_{k+1} - z_k) < steptol.
%   Default: maxit = 50, restol = 1e-14, steptol = 1e-14.
%   If maxit, restol, steptol is empty, the default is used.
%
%   [zk, numiter] = NEWTON(f, df, z0) is the number of iterations.
%   Points that have not converged have numiter = maxit + 1.

% Set default values:
if ((nargin < 4) || isempty(maxit)), maxit = 50; end
if ((nargin < 5) || isempty(restol)), restol = 1e-14; end
if ((nargin < 6) || isempty(steptol)), steptol = 1e-14; end
% Setup
active = 1:numel(z);                    % track which points to iterate
numiter = (maxit+1)*ones(size(z));      % number of iterations
% Check initial points:
fz = f(z);
converged = (abs(fz) < restol);         % already converged
diverged = (isnan(fz) | isinf(fz));     % NaN or Inf
numiter(converged) = 0;                 % no Newton steps required
active(converged | diverged) = [];      % do not further iterate
fz = fz(active);
% Harmonic Newton iteration:
for kk = 1:maxit
    if (isempty(active)), return, end   % every point has converged
    zold = z(active);
    % Newton step:
    z(active) = zold - df(zold).\fz;
    % Convergence check:
    fz = f(z(active));
    converged = (abs(fz) < restol) | ...
        (abs(z(active) - zold) < steptol*abs(zold));
    diverged = (isnan(fz) | isinf(fz)); % NaN or inf
    numiter(active(converged)) = kk;    % took kk iterations
    active(converged | diverged) = [];  % remove converged points
    fz(converged | diverged) = [];      % remove converged points
end
end
