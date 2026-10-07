function m = measurement_inputs_config()
m=measurement_robustness_config();
m.case_ids={'domestic_inputs'};
m.make_plots=false; 
d=measurement_parameter_data();
k=find(strcmp({m.cases.id},'domestic_inputs'));
assert(numel(k)==1,'measurement:Case','Domestic-input case missing');
m.cases(k).shares=d.domestic_inputs;
m.cases(k).targets=[];
end
