function [tpath, tpath_pts, crossings] = transport_path(fun, crit, caus, segment, scaling)
%TRANSPORT_PATH   Transport path for the transport of images method.
%   tpath = transport_path(fun, crit, caus, [a, b]) with function struct
%   fun, critical set crit and caustics caus, constructs a transport path
%   from the point a to b.  tpath is a cell array of transport steps as
%   described in tstep.
% 
%   tpath = transport_path(fun, crit, caus, [a, b], scaling) prescribes the
%   factor 0 < scaling < 1/2, by which the points for a caustic crossing
%   are constructed.  (Default: scaling = 1/4.)
% 
%   [tpath, tpath_pts, crossings] = transport_path(...) also returns the
%   list of points on the transport path, and a matrix with information on
%   the caustic crossings from caustic_intersection.
%  
% See also tstep, caustic_intersection.

% Set optional scaling parameter:
if ( ( nargin < 5 ) || isempty(scaling) )
    scaling = 1/4;
end

% Get caustic crossings:
crossings = caustic_intersection(fun, crit, caus, segment);
numCausCrossings = size(crossings, 1);

% Distance between crossings:
mindistcross = mindist(crossings(:, 2));
if ( mindistcross < 1e-8 )
    warning("transport_path:crossing", "Minimal distance between crossings is " + mindistcross)
end
% Absolute value of tangent (closeness to cusps):
z0 = crossings(:, 1);
dhdg = fun.dh(z0) .* fun.dg(z0);
abspsi = 2 * abs(imag(sqrt(dhdg).*dhdg ./ ...
    (fun.ddg(z0) .* fun.dh(z0) - fun.dg(z0) .* fun.ddh(z0))));
if ( min(abspsi) < 1e-10 )
    warning("transport_path:cusp", "Crossing close to cusp? abs(tangent) = " + min(abspsi) + ".")
end

% Endpoints of the segment:
a = segment(1);
b = segment(2);

% Construct transport path:
if ( numCausCrossings == 0 )
    tpath{1} = tstep(a, b,0);
    tpath_pts = [a; b];
else
    % Relevant points to construct the transport path:
    landmark = [a; crossings(:,2); b];
    
    % Minimal distance from each crossing point to next landmark points:
    distance = zeros(numCausCrossings, 1);
    for jj = 1:numCausCrossings
        distance(jj) = min(abs(landmark(jj+1) - landmark([jj; jj+2])));
    end
    % Scaling:
    distance = scaling * distance;
    
    % Points on transport path:
    tpath_pts = zeros(2*numCausCrossings + 2, 1);
    tpath_pts(1) = a;
    for jj = 1:numCausCrossings
        tpath_pts([2*jj, 2*jj+1]) = landmark(jj+1) + distance(jj) * sign(b-a) * [-1; 1];
        %if crossings(jj,3) > 0
        %    tpath_pts([2*jj, 2*jj+1]) = landmark(jj+1) + distance(jj) * sign(crossings(jj,4)) * [-1; 1];
        %else
        %    tpath_pts([2*jj, 2*jj+1]) = landmark(jj+1) + distance(jj) * sign(crossings(jj,4)) * [1; -1];
        %end
    end
    tpath_pts(end) = b;
    
    % List of transport steps:
    tpath = cell(2*numCausCrossings + 1, 1);
    for jj = 1:numCausCrossings
        tpath{2*jj-1} = tstep(tpath_pts(2*jj-1), tpath_pts(2*jj),0);
        tpath{2*jj} = tstep(tpath_pts(2*jj), tpath_pts(2*jj+1),0, ...
            crossings(jj, :));
    end
    tpath{end} = tstep(tpath_pts(end-1), tpath_pts(end),0);
end

end