function cfg = solver_config()
root=fileparts(fileparts(mfilename('fullpath')));

cfg.model='paper'; 
cfg.horizons=200; 
cfg.validation_horizons=200; 
cfg.robustness_horizons=[50 75 100 200]; 

cfg.p=struct();
cfg.p.beta=0.96;
cfg.p.delta=0.05;
cfg.p.epsilon=0.22;
cfg.p.xi=1.0;
cfg.p.omega_m=0.1;
cfg.p.omega_s=0.9;
cfg.p.zeta_m=0.53288122606394106;
cfg.p.zeta_s=-cfg.p.zeta_m; 

cfg.p.alpha_m=0.4146;
cfg.p.beta_m=0.3024;
cfg.p.alpha_s=0.4452;
cfg.p.beta_s=0.1035;
cfg.p.alpha_x=0.2223;
cfg.p.beta_x=0.5415;
cfg.p.alpha_e=0.3506;
cfg.p.beta_e=0.4737;

cfg.p.g_A_m = 1.033128;
cfg.p.g_A_s = 0.995785;
cfg.p.g_A_x = 1.072178;
cfg.p.g_A_e = 1.011804;
cfg.A0=[0.9999999999999989; 1.0; 1.0000000000000002; 1.0000000000000002];

cfg.targets=struct('K0',90.457665125476964,'R0',778219.55318695167);

cfg.policies=policy_config(cfg.model);
cfg.policy_ids={'service_temporary','service_permanent', ...
    'consumption_temporary','consumption_permanent','carbon_temporary', ...
    'investment_purchase_20','investment_purchase_matched'};
cfg.policies=cfg.policies(ismember({cfg.policies.id},cfg.policy_ids));

cfg.warm_start_file=fullfile(root,'data','paper_initial_paths.mat');

cfg.scale_h0=5.9911583164995754e-05;
cfg.scale_E0=102.59573262545696;

cfg.newton_options=struct(); 
cfg.tolerance=1e-9;
cfg.make_plots=false;
cfg.output_root=fullfile(root,'outputs');
cfg.run_kind='experiments'; 
cfg.output_tag=''; 
cfg.audit.derivatives=false;
cfg.audit.restart=false;
cfg.audit.horizons=false;
end
