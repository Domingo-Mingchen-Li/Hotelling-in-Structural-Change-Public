function same = same_run_metadata(a,b)
if isequaln(a,b),same=true;return;end
same=false;
if isstruct(a)&&isstruct(b)&&isequal(size(a),size(b))
    names=sort(fieldnames(a));if ~isequal(names,sort(fieldnames(b))),return;end
    for i=1:numel(a)
        for j=1:numel(names)
            if ~same_run_metadata(a(i).(names{j}),b(i).(names{j})),return;end
        end
    end
    same=true;
elseif isnumeric(a)&&isnumeric(b)&&isequal(size(a),size(b))
    da=double(a);db=double(b);scale=max(1,max(abs(da),abs(db)));
    same=all(isfinite(da(:)))&&all(isfinite(db(:))) ...
        &&all(abs(da(:)-db(:))<=32*eps(scale(:)));
end
end
