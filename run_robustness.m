setup_solver;
cfg=solver_config;
cfg.run_kind='robustness';
cfg.horizons=cfg.robustness_horizons;
cfg.audit.horizons=true;
cfg.audit.derivatives=true;
cfg.audit.restart=false;
report=run_solver(cfg);
