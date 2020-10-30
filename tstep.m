function step = tstep(eta0, eta1, depth, varargin)
%TSTEP   Constructs a transport step object/structure.
%   step = TSTEP(eta0, eta1, depth) construct a transport step from eta0 to
%   eta1 on recursion depth 'depth'.
% 
%   step = TSTEP(eta0, eta1, depth, crossingInfo) constructs a transport
%   step where a caustic crossing occurs.  crossingInfo is a row vector
%   with the information from the crossings matrix from
%   caustic_intersection.
% 
% See also caustic_intersection.

step.initial = eta0;
step.end = eta1;
step.depth = depth;
step.crossingInfo = [];

if ( nargin == 4 )
    step.crossingInfo = varargin{1};
end

end
