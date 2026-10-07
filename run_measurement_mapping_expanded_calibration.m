function reports = run_measurement_mapping_expanded_calibration(case_id,resume_file)
setup_solver();m=measurement_mapping_expanded_config();
if nargin<1,case_id='';end
if nargin<2,resume_file='';end
reports=run_measurement_robustness_calibration(case_id,resume_file,m);
end
