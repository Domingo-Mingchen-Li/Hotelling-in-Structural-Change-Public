function c = calibration_config()
root=fileparts(fileparts(mfilename('fullpath')));
c.base=solver_config();
c.base.p.g_A_m=1.033128;
c.base.p.g_A_s=0.995785;
c.base.p.g_A_x=1.072178;
c.base.p.g_A_e=1.011804;
c.base.warm_start_file=fullfile(root,'data','paper_initial_paths.mat');
c.T=200;
c.check_T=300; 
c.moment_names={'capital_output_t0','service_labour_t0','service_labour_change_0_14'};
c.targets=[1.63;0.775083520605178;0.07199915550266978];
c.scales=[1.63;1;0.10]; 
c.fit_tolerances=[1e-5;1e-6;1e-6]; 
c.starts=[c.base.p.zeta_m,c.base.targets.K0,c.base.targets.R0;
 0.9*c.base.p.zeta_m,0.9*c.base.targets.K0,1.1*c.base.targets.R0;
 1.1*c.base.p.zeta_m,1.1*c.base.targets.K0,0.9*c.base.targets.R0];
c.lower=[0,1,1e3]; c.upper=[5,1e5,1e10];
c.start_agreement_tolerance=1e-3; 
c.max_iterations=60;
c.fd_step=1e-3; 
c.max_step=0.25; c.minimum_step=1e-7;
c.initial_damping=1e-3; c.max_damping=1e12;
c.continuation_initial_step=0.25; c.continuation_minimum_step=1/1024;
c.jacobian_relative_tolerance=0.05;
c.maximum_condition_number=1e6;
c.horizon_tolerances=[1e-3;1e-4;1e-4];
c.output_root=c.base.output_root;
end
