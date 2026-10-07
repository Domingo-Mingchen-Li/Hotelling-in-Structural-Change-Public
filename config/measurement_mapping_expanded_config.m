function m = measurement_mapping_expanded_config()
m=measurement_mapping_config();
root=fileparts(fileparts(mfilename('fullpath')));
m.output_root=fullfile(root,'outputs','measurement_robustness_mapping_expanded');
m.case_ids={'mapping_c26'}; 
m.calibration_seed_file=fullfile(root,'data','measurement', ...
 'mapping_c26_checkpoint_seed.mat');
raw=load(m.calibration_seed_file,'warm_starts');
w=select_initial_path(raw.warm_starts,'paper','baseline',m.T);
assert(w.passed,'measurement:Seed','Checkpoint equilibrium did not pass');
p=[w.p.zeta_m,w.Targets.K0,w.Targets.R0];
m.calibration.lower=[0,0.01,1];
m.calibration.starts=[p; p.*[0.9,0.9,1.1]; p.*[1.1,1.1,0.9]];
m.calibration.display_inner=false;
end
