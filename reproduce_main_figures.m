function reproduce_main_figures()
setup_solver();c=paper_figure_config();
export_paper_figures(c.transition,'',c.output,'main');
export_counterfactual_writing(c.results,c.cf);
export_policy_evaluation_writing(c.results,c.welfare);
export_paper_cumulative_table(c.transition,c.output);
fprintf('Figures1-16: %s\n',c.output);
end
