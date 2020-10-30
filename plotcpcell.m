function varargout = plotcpcell(ptscell, varargin)
%PLOTCPCELL   plotcp for each entry of a cell array of points.

h = cell(numel(ptscell), 1);

if ( ishold )
    for jj = 1:numel(ptscell)
        h{jj} = plotcp(ptscell{jj}, varargin{:});
    end
else
    h{1} = plotcp(ptscell{1}, varargin{:});
    hold on
    for jj = 2:numel(ptscell)
        h{jj} = plotcp(ptscell{jj}, varargin{:});
    end
    hold off
end

if ( nargout == 1 )
    varargout{:} = h;
end

end