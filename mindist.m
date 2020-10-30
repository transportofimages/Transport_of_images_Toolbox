function dist = mindist(z)
%MINDIST   Minimal distance of points.
%   dist = mindist(z) is the minimal distance between the entries of the
%   vector z.

if ( numel(z) == 1 )
    % inf { |z_i - z_j| } = inf emptyset = inf.
    dist = inf;
    return
end

% Matrix of distances:
M = abs(bsxfun(@minus, z, z.'));

% The diagonal is zero, replace by NaN, which is ignored by min.
M = M + diag(nan(length(z), 1)); 

% Find min:
dist = min(M, [], "all");

end
