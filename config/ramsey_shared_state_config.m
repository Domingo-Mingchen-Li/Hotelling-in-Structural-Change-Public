function w=ramsey_shared_state_config()
root=fileparts(fileparts(mfilename('fullpath')));
w.reference_report_file=''; 
w.reference_output_root=fullfile(root,'outputs','ramsey');
w.climate_report_file=''; 
w.output_root=fullfile(root,'outputs','ramsey_shared_state_counterfactual');
w.common_beta=[]; 
w.production_initial_step=0.10;w.production_maximum_step=0.25;
w.production_minimum_step=1/1024;w.step_growth=1.5;
w.maximum_continuation_attempts=100;
w.display_production_inner=false;
w.ramsey=ramsey_config();w.ramsey.tools={'resource'};
w.ramsey.check_horizon=true;
w.ramsey.inner_options.residual_tolerance=1e-12; 
w.make_comparison_plots=false;
end
