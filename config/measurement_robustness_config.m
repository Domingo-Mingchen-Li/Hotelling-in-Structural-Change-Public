function m = measurement_robustness_config()
root=fileparts(fileparts(mfilename('fullpath')));
m.base=solver_config();
assert(strcmp(m.base.model,'paper'),'measurement:Mode','Measurement workflow requires paper mode');
m.seed_file=fullfile(root,'data','paper_initial_paths.mat');
m.base.warm_start_file=m.seed_file;
m.output_root=fullfile(root,'outputs','measurement_robustness');
m.T=200;
m.make_plots=false;
m.audit_derivatives=true;
m.audit_restart=true;
m.production_continuation_step=0.25;
m.production_continuation_minimum_step=1/1024;
m.case_ids={'annual_average'};
m.policy_ids=m.base.policy_ids;
m.cumulative_periods=[10,50,100,200]; 

d=measurement_parameter_data();
m.cases(1)=struct('id','annual_average','description', ...
 'Mean of complete annual alpha/beta pairs, global inputs, 2000-2014', ...
 'shares',d.annual_average,'targets',[], ...
 'provenance','data/measurement/measurement_parameters.json');
m.cases(2)=struct('id','domestic_inputs','description', ...
 'Domestic inputs; original pooled CAP/(CAP+COMP); mean annual beta', ...
 'shares',d.domestic_inputs,'targets',[], ...
 'provenance','data/measurement/measurement_parameters.json');

m.calibration=calibration_config();
m.calibration.base=m.base; 
m.calibration.T=m.T;
m.calibration.check_T=300; 
end
