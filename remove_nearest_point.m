function [list, mindist, pos] = remove_nearest_point(list, point)
%REMOVE_NEAREST_POINT   Remove a closest point from a vector.
%   list = remove_nearest_point(list, point) removes the entry from list,
%   which is closest to point.
% 
%   [list, mindist, pos] = remove_nearest_point(list, point) also returns
%   the minimal distance between point and list, and the index pos in list
%   of a nearest element.

% Find closest points:
[mindist, I] = min(abs(list - point));

% Remove one closest point:
list(I(1)) = [];

if ( nargout == 3 )
    pos = I(1);
end

end
