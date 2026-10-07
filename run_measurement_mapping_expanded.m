function reports = run_measurement_mapping_expanded(case_id)
setup_solver();m=measurement_mapping_expanded_config();
if nargin<1,case_id='';end
reports=run_measurement_robustness(case_id,m);
end
