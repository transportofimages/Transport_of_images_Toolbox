% Example from README.md of the Transport of images toolbox

rng(1);
fun = harmonicRat([1 0 0], 1, [2 1], [1 1 0], 2, 0);
zer = tiroots(fun)
res = max(abs(fun.f(zer)))