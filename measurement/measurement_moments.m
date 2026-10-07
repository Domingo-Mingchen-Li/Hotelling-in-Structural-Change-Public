function v = measurement_moments(path)
assert(numel(path.K)>=15,'measurement:Moments','Need dates 0:14');
share=path.l_s([1,15])./(path.l_m([1,15])+path.l_s([1,15]));
v=[path.K(1)/(path.E(1)+path.I(1));share(1);share(2)-share(1)];
assert(all(isfinite(v)),'measurement:Moments','Invalid baseline moments');
end
