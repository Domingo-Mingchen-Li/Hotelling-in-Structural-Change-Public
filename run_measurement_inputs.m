function report = run_measurement_inputs()
setup_solver();m=measurement_inputs_config();
report=run_measurement_robustness('domestic_inputs',m);
end
