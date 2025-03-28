%% (Pflueger) Legacy IRF Plotting function, replicating Figure 3 of Pflueger (2022)
% Plot: 2x3 panel of impulse responses to Short-Term MP shock.
% Subplot 1: MP Rate
% Subplot 2: Output Gap
% Subplot 3: Inflation
% Subplot 4: Equity Return (PD)
% Subplot 5: (Nominal) Bond Yield
% Subplot 6: (Real) Bond Yield

function plot_StructuralIRF_STMP(asset_p, figNameFile, figNameHeader)
        rng(0);
        asset = asset_p;

        % Plot figure
        h = figure;
        pause(0.00001);

        set(gcf,'WindowState','Maximized');
        set(gcf,'color','w');
        x_axis = [0:19];

        subplot(2,3,1);
        p = plot(x_axis,asset.Irf3.i(1:numel(x_axis)), 'linewidth',2);
        set(p,'Color', [0,0,0]);
        hold on;
        xlim([0,numel(x_axis)-1]);
        plot(x_axis,0*asset.Irf3.i(1:numel(x_axis)),'-k');
        set(gca, 'FontSize', 20);
        ylabel('MP Rate','fontweight','normal','fontsize',20);
        
        subplot(2,3,2);
        p = plot(x_axis,asset.Irf3.x(1:numel(x_axis)), 'linewidth',2);
        set(p,'Color', [0,0,0]);
        hold on;
        xlim([0,numel(x_axis)-1]);
        plot(x_axis,0*asset.Irf3.x(1:numel(x_axis)),'-k');
        set(gca, 'FontSize', 20);
        ylabel('Output Gap','fontweight','normal','fontsize',20);
        
        subplot(2,3,3);
        p = plot(x_axis,asset.Irf3.pi(1:numel(x_axis)), 'linewidth',2);
        set(p,'Color', [0,0,0]);
        hold on;
        xlim([0,numel(x_axis)-1]);
        plot(x_axis,0*asset.Irf3.pi(1:numel(x_axis)),'-k');
        set(gca, 'FontSize', 20);
        ylabel('Inflation','fontweight','normal','fontsize',20);

        subplot(2,3,4);
        set(gcf,'color','w');
        p = plot(x_axis,asset.Irf3.PD(1:numel(x_axis)), 'linewidth',2);
        set(p,'Color', [0,0,0]);
        hold on;
        p = plot(x_axis,asset.Irf3.PD_rn(1:numel(x_axis)),'--', 'linewidth',2);
        set(p,'Color','red'); 
        hold on;
        p = plot(x_axis,asset.Irf3.PD_rp(1:numel(x_axis)),':', 'linewidth', 2);
        set(p,'Color', 'blue');
        hold on;
        
        xlim([0,numel(x_axis)-1]);
        plot(x_axis,0*asset.Irf3.PD(1:numel(x_axis)),'-k');
        set(gca, 'FontSize', 20);
        ylabel('Equity Return','fontweight','normal','fontsize',20);
  
        subplot(2,3,5);
        p = plot(x_axis,asset.Irf3.y10nom(1:numel(x_axis)), 'linewidth',2);
        set(p,'Color', [0,0,0]);
        hold on;
        p = plot(x_axis,asset.Irf3.y10nom_rn(1:numel(x_axis)),'--', 'linewidth',2);
        set(p,'Color','red');    
        hold on;
        p = plot(x_axis,asset.Irf3.y10nom_rp(1:numel(x_axis)),':', 'linewidth', 2);
        set(p,'Color', 'blue');
        hold on;
        
        plot(x_axis,0*asset.Irf3.y10nom(1:numel(x_axis)),'-k');
        set(gca, 'FontSize', 20); % to increase ticks font size on axes
        xlim([0,numel(x_axis)-1]);
        ylim([-0.06,0.21]);
        ylabel('Nominal Bond Yield','fontweight','normal','fontsize',20);
        
        subplot(2,3,6);
        p = plot(x_axis,asset.Irf3.y10real(1:numel(x_axis)), 'linewidth',2);
        set(p,'Color', [0,0,0]);
        hold on;
        p = plot(x_axis,asset.Irf3.y10real_rn(1:numel(x_axis)),'--', 'linewidth',2);
        set(p,'Color','red');    
        hold on;
        p = plot(x_axis,asset.Irf3.y10real_rp(1:numel(x_axis)),':', 'linewidth', 2);
        set(p,'Color', 'blue');
        hold on;
        
        plot(x_axis,0*asset.Irf3.y10real(1:numel(x_axis)),'-k');

        hold off
        set(gca, 'FontSize', 20); % to increase ticks font size on axes
        xlim([0,numel(x_axis)-1]);
        ylim([-0.06,0.21]);
        ylabel('Real Bond Yield','fontweight','normal','fontsize',20);
        legend({'Overall','Risk Neutral', 'Risk Premium'}, 'Location', 'northeast');
        legend('boxoff');
        
        pause(.3);
        
        saveas(gcf,'IRF_STMP','jpg')
        
end
