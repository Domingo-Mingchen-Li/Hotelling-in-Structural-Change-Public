function plot_Figure1_Baseline()
    clc; clear; close all;

% Input the path to any MATLAB data file (.mat) located in the Sim_Results directory.
    file_path = 'Sim_Results\service_subsidy_10_periods.mat';
    
    obs_period = 40; 
    
    if ~exist(file_path, 'file')
        error('File not found. Please specify the correct data path.');
    end
    load(file_path, 'SimData');
    
    Base = SimData.BaseRes;
    time = 0:obs_period;
    
    K       = Base.K(1:obs_period+1);
    R       = Base.R(1:obs_period+1);
    e       = Base.R_input(1:obs_period+1);
    Pm_Ps   = Base.p_m(1:obs_period+1) ./ Base.p_s(1:obs_period+1);
    
    L_total = Base.l_m(1:obs_period+1) + Base.l_s(1:obs_period+1);
    Lm_L    = Base.l_m(1:obs_period+1) ./ L_total;
    Ls_L    = Base.l_s(1:obs_period+1) ./ L_total;
    
    e_Y     = Base.R_input(1:obs_period+1) ./ Base.Y_GDP(1:obs_period+1);

    figure('Name', 'Figure 1: Baseline Dynamics', 'Color', 'w', 'Position', [100, 100, 1600, 800]);
    
    lw = 2.5; 
    fs = 14;  

    subplot(2, 4, 1);
    plot(time, K, 'b-', 'LineWidth', lw);
    title('Capital Stock (K)', 'FontSize', fs);
    apply_format(obs_period);

    subplot(2, 4, 2);
    plot(time, Pm_Ps, 'm-', 'LineWidth', lw);
    title('Relative Price (Pm / Ps)', 'FontSize', fs);
    apply_format(obs_period);

    subplot(2, 4, 3);
    plot(time, Lm_L, 'b-', 'LineWidth', lw);
    title('Labor Share: Manuf (L_m/L)', 'FontSize', fs);
    apply_format(obs_period); 
    ylim([0, 1]);

    subplot(2, 4, 4);
    plot(time, Ls_L, 'r-', 'LineWidth', lw);
    title('Labor Share: Services (L_s/L)', 'FontSize', fs);
    apply_format(obs_period); 
    ylim([0, 1]);

    subplot(2, 4, 5);
    plot(time, R, 'r-', 'LineWidth', lw);
    title('Resource Stock (R)', 'FontSize', fs);
    apply_format(obs_period);

    subplot(2, 4, 6);
    plot(time, e, 'm-', 'LineWidth', lw);
    title('Extraction Flow (e)', 'FontSize', fs);
    apply_format(obs_period);

    subplot(2, 4, 7);
    plot(time, e_Y, 'm-', 'LineWidth', lw);
    title('Resources Intensity (e / Y)', 'FontSize', fs);
    apply_format(obs_period);

    fprintf('Figure 1 rendered successfully.\n');
end

function apply_format(obs_p)
    set(gca, 'FontSize', 11, 'TickDir', 'out');
    grid on;
    xlim([0, obs_p]);
end