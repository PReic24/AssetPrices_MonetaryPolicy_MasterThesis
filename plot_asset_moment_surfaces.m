%% (Reichert) Plot asset moment surfaces given the monetary policy grid of
% Equity Premium
% Equity Volatility
% Equity Sharpe Ratio
% Nominal Bond log Yield Spread
% Nominal Bond Excess Return Volatility
% 1-Year Excess Returns on log Yield Spread
function plot_asset_moment_surfaces(simulation_results, gamma_x_list, gamma_pi_list)
    simulation_results_plot = simulation_results;

    x = gamma_x_list;
    y = gamma_pi_list;

    moment_names = ["Equity Premium", "Equity Volatility", "Equity Sharpe Ratio", "Yield Spread", "Excess Bond Return Volatility", "1-Year-Return on Yield Spread"];
    
    moments_to_plot.eq_premium                  = simulation_results_plot.equity.eq_premium;
    moments_to_plot.eq_volatility               = simulation_results_plot.equity.vol;
    moments_to_plot.eq_sharpe_ratio             = simulation_results_plot.equity.sharpeRatio;

    moments_to_plot.yield_spread                = simulation_results_plot.nominal_bonds.mean_log_yield_spread;
    moments_to_plot.vol_bond_excess             = simulation_results_plot.nominal_bonds.vol;
    moments_to_plot.yr_yieldspread_coeff        = simulation_results_plot.nominal_bonds.coeffRegRetOnYS1y;

    fields = fieldnames(moments_to_plot);
    
    %% Start plotting the figure
    figure;
    set(gcf, 'WindowState', 'Maximized', 'color', 'w');

    for j = 1:3
        subplot(1, 3, j)
        p = surf(x,y,moments_to_plot.(fields{j}));
        hold on

        % formatting
        title(moment_names{j})
        xlabel("\gamma_{x}", "FontSize", 10)
        ylabel("\gamma_{\pi}", "FontSize", 10)
        hold off

        if j == 1 || j ==2
            view([135,20])
        else
            view([-45, 20])
        end
    end
    exportgraphics(gcf, 'Equity_surfaces.jpg', 'ContentType', 'vector', 'Resolution', 600)

    figure
    set(gcf, 'WindowState', 'Maximized', 'color', 'w');
    
    for k = 1:3
        l = k + 3;
        subplot(1, 3, k)
        p = surf(x,y,moments_to_plot.(fields{l}));
        hold on

        % formatting
        title(moment_names{l})
        xlabel("\gamma_{x}", "FontSize", 10)
        ylabel("\gamma_{\pi}", "FontSize", 10)
        hold off
        
        if l == 4 || l ==5
            view([135,20])
        else
            view([-45, 20])
        end

    end

    exportgraphics(gcf, 'Bond_surfaces.jpg', 'ContentType', 'vector', 'Resolution', 600)
end