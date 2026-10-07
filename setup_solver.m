function root=setup_solver()
root=fileparts(mfilename('fullpath'));
addpath(root);
folders={'core','config','run','audit','analysis','calibration','measurement', ...
 'data/measurement','counterfactual','climate','ramsey','ramsey_shared'};
for k=1:numel(folders),addpath(fullfile(root,folders{k}));end
end
