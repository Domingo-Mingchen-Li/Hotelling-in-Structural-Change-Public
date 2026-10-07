function m = climate_config()
root=fileparts(fileparts(mfilename('fullpath')));
m.base=solver_config();
m.T=200; m.check_T=300;
m.L=10;                         
m.reference_date=50;             
m.initial_components=zeros(4,1); 
m.a=[0.2173;0.2240;0.2824;0.2763];
m.timescales_years=[Inf;394.4;36.54;4.304];
m.kernel_source='Joos et al. (2013), ACP 13:2793-2825, Table 5';
m.kernel_source_url='https://doi.org/10.5194/acp-13-2793-2013';
m.kernel_fit_years=1000;
m.etas=[2 1]; m.loss_targets=[0.005 0.01 0.02];
m.benchmark_eta=2; m.benchmark_loss=0.01;
m.loss_definition='proportional_basket';
m.maximum_utility_tail_years=3000;
m.utility_tail_block_years=100;
m.utility_tail_relative_tolerance=1e-10;
m.horizon_relative_tolerance=5e-4;
m.identity_tolerance=1e-9;
m.baseline_tail_gap_tolerance=1e-3;
m.display_inner=false;
m.audit_derivatives=true; m.audit_restart=false;
m.seed_file=fullfile(root,'data','paper_initial_paths.mat');
m.output_root=fullfile(root,'outputs','climate_parameterization');
end
