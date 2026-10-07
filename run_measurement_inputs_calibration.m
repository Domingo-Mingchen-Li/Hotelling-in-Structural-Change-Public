function report = run_measurement_inputs_calibration(resume_file)
setup_solver();m=measurement_inputs_config();
if nargin<1,resume_file='';end
report=run_measurement_robustness_calibration('domestic_inputs',resume_file,m);
end
