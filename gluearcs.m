function curves = gluearcs(arcs, tol)
%GLUEARCS   Glue critical arcs together.
%   CC = GLUEARCS(ARCS) takes the n by k matrix ARCS containing k critical
%   arcs and glues consecutive arcs together.  CC is a cell array, where
%   each cell contains one critical curve.
%
%   CC = GLUEARCS(ARCS, tol) specifies the tolerance.  Default tol = 1e-5.

if ( ( nargin < 2 ) || isempty(tol) )
    tol = 1e-5;
end

% Begin first critical curve:
curves{1} = arcs(:,1);
endpoint = arcs(end,1);
activeCols = 2:size(arcs, 2);

while ( ~isempty(activeCols) )
    [mm, I] = min(abs(endpoint - arcs(1,activeCols)));
    if ( mm < tol )
        % Add arc to current curve: By construction in critical_curves,
        % first entry of new column = last entry of old column, hence glue
        % only the entries 2:end.
        curves{end} = [curves{end}; arcs(2:end,activeCols(I))];
        endpoint = arcs(end, activeCols(I));
        activeCols(I) = [];
    else
        % Begin new critical curve:
        curves{end+1} = arcs(:,activeCols(1));
        endpoint = arcs(end, activeCols(1));
        activeCols(1) = [];
    end
end

end
