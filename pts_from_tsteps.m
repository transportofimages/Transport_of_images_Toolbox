function pts = pts_from_tsteps(steps)
%PTS_FROM_TSTEPS   List of points from list of transport steps.
%   pts = pts_from_tsteps(steps) are the points eta_k from a list of
%   transport steps.

if ( isempty(steps) )
    pts = [];
    return
end

pts = zeros(length(steps) + 1, 1);
pts(1) = steps{1}.initial;
for jj = 1:length(steps)
    pts(jj + 1) = steps{jj}.end;
end

end