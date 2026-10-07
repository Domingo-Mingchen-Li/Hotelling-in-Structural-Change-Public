setup_solver;
cfg=solver_config;
cfg.run_kind='validation';
cfg.horizons=cfg.validation_horizons;
cfg.audit.derivatives=true;
cfg.audit.restart=true;
report=run_solver(cfg);
