function [sol, sol_hist, prediction_hist, tpath_pts, numsol, total_steps, numiter_newton, tsteps_list] ...
    = transport_phase_hist(fun, sol, tpath, tpath_pts, ...
    tol_newton_accuracy, tol_sep, maxit_newton, max_depth)
%TRANSPORT_PHASE_HIST   Transport phase.
%   sol = transport_phase_hist(fun, sol, tpath) transports the solutions of
%   fun.f = eta_k for eta_k on the transport path tpath.  If the transport
%   is not successfull, sol = NaN is returned.
% 
%   sol = transport_phase_hist(fun, sol, tpath, tol_newton_accuracy, ...
%   tol_sep, maxit_newton, maxnum_tsteps) allows to specify parameters:
%   tol_newton_accuracy is the accuracy when evaluating residuals in the
%   harmonic Newton method (default 1e-10), tol_sep is the minimal distance
%   between two solutions to be considered distinct (default 1e-10),
%   maxit_newton is the maximal number of Newton steps (default 100),
%   max_depth is the maximal recursio depth for the refinement of a
%   transport path.
% 
%   [sol, sol_hist, prediction_hist] = transport_phase_hist(...) collects the solutions in
%   each step of the transport phase in the cell array sol_hist.

% Set defauls if not provided:
if ( ( nargin < 5 ) || isempty(tol_newton_accuracy) )
    tol_newton_accuracy = 1e-10;
end
if ( ( nargin < 6 ) || isempty(tol_sep) )
    tol_sep = 1e-10;
end
if ( ( nargin < 7 ) || isempty(maxit_newton) )
    maxit_newton = 100;
end
if ( ( nargin < 8 ) || isempty(max_depth) )
    max_depth = 10;
end


sol_hist{1} = sol;
prediction_hist = {};

jobs = tpath;
tsteps_list = {};
numsol = length(sol);

% Initialize counters:
total_steps = 0;
numiter_newton = 0;

while ( ~isempty(jobs) )
    % Retrieve next step:
    currentstep = jobs{1};
    jobs(1) = [];
    
    % Count total number of steps.
    total_steps = total_steps + 1;
    % Check recursion depth.
    if ( currentstep.depth > max_depth )
        sol = NaN;
        return;
    end
    
    % Perform the homotopy step.
    if ( isempty(currentstep.crossingInfo) )
        % Step without caustic crossing.
        
        % Newton:
        [newsol, numiter] = harmonicNewton(@(z) fun.f(z) - currentstep.end, ...
            fun.dh, fun.dg, sol, maxit_newton);
        numiter_newton = numiter_newton + sum(numiter, "all");
        
        % Check if the new solutions are accurate and distinct:
        if ( ( norm(fun.f(newsol) - currentstep.end, inf) < tol_newton_accuracy ) && ...
                ( mindist(newsol) >= tol_sep ) )
            % Transport step was successful:
            prediction_hist{end+1,1} = sol;
            sol = newsol;
            % numsol stays the same.
            sol_hist{end+1,1} = sol;
            tsteps_list{end+1, 1} = currentstep;
        else
            % Divide step in two:
            midpt = (currentstep.initial + currentstep.end)/2;
            jobs = [{tstep(currentstep.initial, midpt,currentstep.depth+1)};
                {tstep(midpt, currentstep.end, currentstep.depth+1)};
                jobs];
            I = find(tpath_pts == currentstep.initial);
            tpath_pts = [tpath_pts(1:I); midpt; tpath_pts(I+1:end)];
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
        sqrt_term = sqrt(t*complex_dilat);
        init_pts = [z0 + 1i*sqrt_term; z0 - 1i*sqrt_term];
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
            error("transport_phase_hist:multiple_crossings", "Multiple crossings not yet implemented.")
        end
        
        % Newton:
        [newsol, numiter] = harmonicNewton(@(z) fun.f(z) - currentstep.end, ...
            fun.dh, fun.dg, prediction, maxit_newton);
        numiter_newton = numiter_newton + sum(numiter, "all");
        
        % Check if the new solutions are accurate and distinct:
        if ( ( norm(fun.f(newsol) - currentstep.end, inf) < tol_newton_accuracy ) && ...
                ( mindist(newsol) >= tol_sep ) )
            % Transport step was successful:
            sol = newsol;
            numsol = numsol + numsolchange;
            sol_hist{end+1,1} = sol;
            prediction_hist{end+1,1} = prediction;
            tsteps_list{end+1, 1} = currentstep;
        else
            % Divide step in three:
            intermed1 = (3*currentstep.initial + currentstep.end)/4;
            intermed2 = (currentstep.initial + 3*currentstep.end)/4;
            jobs = [{tstep(currentstep.initial, intermed1, currentstep.depth+1)};
                {tstep(intermed1, intermed2, currentstep.depth+1, currentstep.crossingInfo)};
                {tstep(intermed2, currentstep.end, currentstep.depth+1)};
                jobs];
            I = find(tpath_pts == currentstep.initial);
            tpath_pts = [tpath_pts(1:I); intermed1; intermed2; tpath_pts(I+1:end)];
        end
    end
end

end