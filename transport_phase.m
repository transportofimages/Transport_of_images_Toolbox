function [sol, numsol, steps_total, numiter_newton, steps_failed, tsteps_list] = transport_phase(fun, sol, tpath, ...
    tol_newton_accuracy, tol_sep, maxit_newton, max_depth)
%TRANSPORT_PHASE   Transport phase.
%   sol = transport_phase(fun, sol, tpath) transports the solutions of
%   fun.f = eta_k for eta_k on the transport path tpath.  If the transport
%   is not successfull, sol = NaN is returned.
% 
%   sol = transport_phase(fun, sol, tpath, tol_newton_accuracy, tol_sep,
%   ... maxit_newton, maxnum_tsteps) allows to specify parameters:
%   tol_newton_accuracy is the accuracy when evaluating residuals in the
%   harmonic Newton method (default 1e-10), tol_sep is the minimal distance
%   between two solutions to be considered distinct (default 1e-10),
%   maxit_newton is the maximal number of Newton steps (default 100),
%   max_depth is the maximal recursio depth for the refinement of a
%   transport path.

% Set defauls if not provided:
if ( ( nargin < 4 ) || isempty(tol_newton_accuracy) )
    tol_newton_accuracy = 1e-10;
end
if ( ( nargin < 5 ) || isempty(tol_sep) )
    tol_sep = 1e-10;
end
if ( ( nargin < 6 ) || isempty(maxit_newton) )
    maxit_newton = 100;
end
if ( ( nargin < 7 ) || isempty(max_depth) )
    max_depth = 10;
end

jobs = tpath;
tsteps_list = {};
numsol = length(sol);

% Initialize counters:
steps_failed = 0;
steps_total = 0;
numiter_newton = 0;

while ( ~isempty(jobs) )
    % Retrieve next step:
    currentstep = jobs{1};
    jobs(1) = [];
    
    % Check recursion depth.
    if ( currentstep.depth > max_depth )
        sol = NaN;
        return;
    end
    
    % Count total number of steps.
    steps_total = steps_total + 1;
    
    % Perform the homotopy step.
    if ( isempty(currentstep.crossingInfo) )
        % Step without caustic crossing.
        
        % Newton:
        [newsol, numiter] = harmonicNewton(@(z) fun.f(z) - currentstep.end, ...
            fun.dh, fun.dg, sol, maxit_newton);
        numiter_newton = numiter_newton + sum(numiter, "all");
        
        % Check if the new solutions are accurate and distinct:
        if ( ( numsol == 0 ) || ...
                ( ( norm(fun.f(newsol) - currentstep.end, inf) < tol_newton_accuracy ) && ...
                ( mindist(newsol) >= tol_sep ) ) )
            % Transport step was successful:
            sol = newsol;
            % numsol stays the same.
            tsteps_list{end+1, 1} = currentstep;
        else
            % Divide step in two:
            midpt = (currentstep.initial + currentstep.end)/2;
            jobs = [{tstep(currentstep.initial, midpt, currentstep.depth + 1)};
                {tstep(midpt, currentstep.end, currentstep.depth + 1)};
                jobs];
            steps_failed = steps_failed + 1;
        end
        
    else
        % Step with caustic crossing:
        z0 = currentstep.crossingInfo(1);
        fz0 = currentstep.crossingInfo(2);
        numsolchange = currentstep.crossingInfo(3);
        c = currentstep.crossingInfo(4);
        complex_dilat = currentstep.crossingInfo(5); % 1st dilatation
        
        % Determine t:  eta_2 = f(z0) + t*c
        t = (currentstep.end - fz0)/c;
        init_pts = z0 + 1i*sqrt(t*complex_dilat) * [1; -1];
        % TODO: should t be made real?
        
        % Construct prediction set:
        if ( numsolchange == 2 )
            % Add two points to the prediction set:
            prediction = [sol; init_pts];
        elseif ( numsolchange == -2 )
            % Remove two points from sol:
            [init_pts, numiter] = harmonicNewton(@(z) fun.f(z) - currentstep.initial, ...
                fun.dh, fun.dg, init_pts, maxit_newton);
            numiter_newton = numiter_newton + sum(numiter, "all");
            prediction = remove_nearest_point(sol, init_pts(1));
            prediction = remove_nearest_point(prediction, init_pts(2));
        else
            error("transport_phase:multiple_crossings", "Multiple crossings not yet implemented.")
        end
        
        % Newton:
        [newsol, numiter] = harmonicNewton(@(z) fun.f(z) - currentstep.end, ...
            fun.dh, fun.dg, prediction, maxit_newton);
        numiter_newton = numiter_newton + sum(numiter, "all");
        
        % Check if the new solutions are accurate and distinct:
        if ( ( numsol + numsolchange == 0 ) || ...
                ( ( norm(fun.f(newsol) - currentstep.end, inf) < tol_newton_accuracy ) && ...
                ( mindist(newsol) >= tol_sep ) ) )
            % Transport step was successful:
            sol = newsol;
            numsol = numsol + numsolchange;
            tsteps_list{end+1, 1} = currentstep;
        else
            % Divide step in three:
            intermed1 = (3*currentstep.initial + currentstep.end)/4;
            intermed2 = (currentstep.initial + 3*currentstep.end)/4;
            jobs = [{tstep(currentstep.initial, intermed1,currentstep.depth + 1)};
                {tstep(intermed1, intermed2,currentstep.depth + 1, currentstep.crossingInfo)};
                {tstep(intermed2, currentstep.end,currentstep.depth + 1)};
                jobs];
            steps_failed = steps_failed + 1;
        end
    end
end

end