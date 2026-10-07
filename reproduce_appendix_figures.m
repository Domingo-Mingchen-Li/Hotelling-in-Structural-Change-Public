function reproduce_appendix_figures()
setup_solver();c=paper_figure_config();
export_paper_figures(c.transition,c.horizon,c.appendix_output,'appendix');
export_measurement_writing(fullfile(c.results,'measurement'),c.measurement);
o=c.welfare;o.destination=c.appendix_output;o.make_plots=false;o.make_tables=true;
export_policy_evaluation_writing(c.results,o);
fprintf('Appendix figures and tables: %s\n',c.appendix_output);
end
