function w=industrial_welfare_config()
root=fileparts(fileparts(mfilename('fullpath')));
w.climate_report_file='';
w.climate_output_root=fullfile(root,'outputs','climate_parameterization');
w.output_root=fullfile(root,'outputs','industrial_welfare');
w.check_horizon=true;             
w.display_inner=false;
w.audit_derivatives=true;w.audit_restart=false;
w.horizon_anchor_tolerance=5e-4;  
w.neutrality_tolerance=1e-7;
w.sign_anchor_tolerance=1e-8;     
w.make_plots=false;
s=struct('id','','start',0,'duration',10,'permanent',false, ...
 'tau_c_m',0,'tau_c_s',0,'tau_e',0,'tau_int',0,'tau_x',0, ...
 'investment_policy_kind','none');
w.policies=repmat(s,1,6);
w.policies(1).id='service_temporary';w.policies(1).tau_c_s=-0.2;
w.policies(2).id='service_permanent';w.policies(2).tau_c_s=-0.2;w.policies(2).permanent=true;
w.policies(3).id='consumption_temporary';w.policies(3).tau_c_m=0.2;w.policies(3).tau_c_s=0.2;
w.policies(4).id='consumption_permanent';w.policies(4).tau_c_m=0.2;w.policies(4).tau_c_s=0.2;w.policies(4).permanent=true;
w.policies(5).id='investment_purchase_matched';w.policies(5).tau_x=-1/6;w.policies(5).investment_policy_kind='purchase_subsidy';
w.policies(6).id='investment_purchase_20';w.policies(6).tau_x=-0.2;w.policies(6).investment_policy_kind='purchase_subsidy';
end
