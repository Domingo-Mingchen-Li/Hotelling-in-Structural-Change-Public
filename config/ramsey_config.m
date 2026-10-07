function w=ramsey_config()
root=fileparts(fileparts(mfilename('fullpath')));
w.climate_report_file=''; 
w.climate_output_root=fullfile(root,'outputs','climate_parameterization');
w.output_root=fullfile(root,'outputs','ramsey');
w.tools={'resource'}; 
w.eta=2;w.loss_target=0.01; 
w.resource_tax_bounds=[0 10];
w.service_wedge_bounds=[-0.3 0.3]; 
w.check_horizon=true; 
w.start_amplitudes=[0 0.1]; 
w.start_decay_years=30;
w.max_iterations=300;w.memory=10;w.max_backtracks=30;
w.maximum_control_step=1;w.armijo=1e-4;
w.projected_gradient_tolerance=2e-5;
w.multistart_gain_tolerance=1e-4; 
w.horizon_gain_tolerance=5e-4;
w.horizon_early_wedge_tolerance=0.005; 
w.horizon_comparison_dates=30;
w.derivative_step=1e-5;
w.gradient_audit_step=1e-4;w.gradient_audit_tolerance=2e-5;
w.audit_gradient=true;w.audit_equilibrium_derivatives=true;
w.display_inner=false; 
w.inner_options=struct('max_iterations',70,'residual_tolerance',1e-10);
w.make_plots=false;
w.terminal_policy='zero_from_T';
end
