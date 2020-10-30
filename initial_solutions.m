function [eta, sol, numiter_newton] = initial_solutions(fun, caus, phi)
%INITIAL_SOLUTIONS   Initial phase of transport of images.
%   [eta, sol] = initial_solutions(fun, caus) computes a sufficiently large
%   initial value eta outside of the caustics caus and all solutions of
%   fun.f(z) = eta.
%
%   [eta, sol] = initial_solutions(fun, caus, phi) determines an initial
%   value eta with angle(eta) = phi. (Default: phi = 0.)
% 
%   [eta, sol, numiter_newton] = initial_solutions(...) returns the overall
%   number of harmonic Newton iterations used.

% Separation tolerance for clustering initial solutions:
tol_sep = 1e-10;

if ( ( nargin < 3 ) || isempty(phi) )
    phi = 0;
end

%% Compute poles, orders and coefficients:
% fun.f(z) = r(z) + conj(s(z))

% fprintf("Compute poles...\n")
polematrix = get_poles(fun);

if ( isempty(polematrix) )
    error("initial_solutions:constant", "Function is constant.")
end

%% Initial solutions

% Initial guess for sufficiently large eta outside of the caustics:
max_on_caus = 0;
for jj = 1:length(caus)
    max_on_caus = max(max_on_caus, max(abs(caus{jj})));
end
eta = 2*max_on_caus * exp(1i*phi);


numiter_newton = 0;

for step = 1:50
    % Find initial solutions:
    % Prediction: construct approximate solutions:
    sol = [];
    for ii = 1:size(polematrix, 1)
        % Approximate solutions at pole ii:
        tmp = approx_zeros_at_pole(polematrix(ii, 1), polematrix(ii, 2), ...
            polematrix(ii, 3), polematrix(ii, 4), polematrix(ii, 5) + eta);
        sol = [sol; tmp];
    end
    
    % Correction:
    [sol, numiter] = harmonicNewton(@(z) fun.f(z) - eta, fun.dh, fun.dg, sol);
    numiter_newton = numiter_newton + sum(numiter, "all");
    
    % Check that solutions are distinct and that the residual is small:
    if ( ( mindist(sol) >= tol_sep ) && ( norm(fun.f(sol) - eta) < 1e-10 ) )
        % Success.
        return
    else
        % Increase eta and try again:
        eta = 2*eta;
    end
end

error("initial_solutions:eta", "Initial phase failed: no suitable eta found.")

end
