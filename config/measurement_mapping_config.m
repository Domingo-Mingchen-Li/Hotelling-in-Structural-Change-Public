function m = measurement_mapping_config()
m=measurement_robustness_config();
m.case_ids={'mapping_c26'};
m.make_plots=false; 
d=sector_mapping_parameter_data();
ids={'mapping_c26'};
descriptions={'C26 moved from investment to manufacturing'};
for i=1:numel(ids)
 id=ids{i};data=d.(id);
 spec=struct('id',id,'description',descriptions{i}, ...
  'shares',data.shares,'targets',data.targets, ...
  'provenance','data/measurement/sector_mapping_parameters.json');
 m.cases(end+1)=spec;
end
end
