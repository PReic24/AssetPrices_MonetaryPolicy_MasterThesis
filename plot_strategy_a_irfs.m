%% (Reichert) Plots impulse response functions given Calibration Strategy A of 
% The Output Gap
% Nominal Interest Rate
% Equity Returns
% 10-Year Nominal Bond Yields
function plot_strategy_a_irfs(simulation_results) 

   
    x_axis = 0:19;

    % Create output figure
    figure;
    set(gcf, 'WindowState', 'Maximized', 'color', 'w');

    base = simulation_results.macro_irfs.x(1:numel(x_axis),1);
    strat_A = simulation_results.macro_irfs.x(1:numel(x_axis),4);
    
    plot(x_axis, base , '-k', 'LineWidth', 2);
    hold on
    plot(x_axis, strat_A ,"--k", 'LineWidth', 2);
    hold on
    plot(x_axis,0*base,'-k')
    hold off
    
    xlim([0, numel(x_axis) - 1]);
    ylabel("Output Gap")
    xlabel("Time Periods")
    set(gca, 'FontSize', 18)

    
    legend('Base', 'Strategy A', 'Location', 'Best');

    set(gca, 'LooseInset', [0,0,0,0]) % Remove excess space
    set(gcf, 'PaperPositionMode', 'auto')

    saveas(gcf, "Output_IRFs_strat_a", 'jpg');

    % Create interest rate figure
    figure;
    set(gcf, 'WindowState', 'Maximized', 'color', 'w');

    base = simulation_results.macro_irfs.i(1:numel(x_axis),1);
    strat_A = simulation_results.macro_irfs.i(1:numel(x_axis),4);
    
    plot(x_axis, base , '-k', 'LineWidth', 2);
    hold on
    plot(x_axis, strat_A ,"--k", 'LineWidth', 2);
    hold on
    plot(x_axis,0*base,'-k')
    hold off
    
    xlim([0, numel(x_axis) - 1]);
    ylabel("Nominal Interest Rate")
    xlabel("Time Periods")
    set(gca, 'FontSize', 18)

    
    legend('Base', 'Strategy A', 'Location', 'Best');

    set(gca, 'LooseInset', [0,0,0,0]) % Remove excess space
    set(gcf, 'PaperPositionMode', 'auto')

    saveas(gcf, "NomInterestRate_IRFs_strat_a", 'jpg');

    % Create equity return figure
    figure;
    set(gcf, 'WindowState', 'Maximized', 'color', 'w');

    base = simulation_results.macro_irfs.pd(1:numel(x_axis),1);
    strat_A = simulation_results.macro_irfs.pd(1:numel(x_axis),4);
    
    plot(x_axis, base , '-k', 'LineWidth', 2);
    hold on
    plot(x_axis, strat_A ,"--k", 'LineWidth', 2);
    hold on
    plot(x_axis,0*base,'-k')
    hold off
    
    xlim([0, numel(x_axis) - 1]);
    ylabel("Equity Return")
    xlabel("Time Periods")
    set(gca, 'FontSize', 18)

    
    legend('Base', 'Strategy A', 'Location', 'Best');

    set(gca, 'LooseInset', [0,0,0,0]) % Remove excess space
    set(gcf, 'PaperPositionMode', 'auto')

    saveas(gcf, "Equity_Return_IRFs_strat_a", 'jpg');

    % Create nominal bond yield figure
    figure;
    set(gcf, 'WindowState', 'Maximized', 'color', 'w');

    base = simulation_results.macro_irfs.y10nom(1:numel(x_axis),1);
    strat_A = simulation_results.macro_irfs.y10nom(1:numel(x_axis),4);
    
    plot(x_axis, base , '-k', 'LineWidth', 2);
    hold on
    plot(x_axis, strat_A ,"--k", 'LineWidth', 2);
    hold on
    plot(x_axis,0*base,'-k')
    hold off
    
    xlim([0, numel(x_axis) - 1]);
    ylabel("Nominal Bond Yield")
    xlabel("Time Periods")
    set(gca, 'FontSize', 18)

    
    legend('Base', 'Strategy A', 'Location', 'Best');

    set(gca, 'LooseInset', [0,0,0,0]) % Remove excess space
    set(gcf, 'PaperPositionMode', 'auto')

    saveas(gcf, "NomBond_Yield_IRFs_strat_a", 'jpg');


    
end
