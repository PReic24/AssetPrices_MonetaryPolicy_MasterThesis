%-------------------------------------------------------------------------%
% Simulation Framework for the master thesis
% Title: On the Relation between Asset Prices and the Stance of Monetary Policy
% Autor: Pascal Reichert
% Names in parenthesis indicate authorship of script / function / code
% structure
%-------------------------------------------------------------------------% 
% Largely based on:
% Minimum working example for programming package
% Autors: Carolin Pﬂueger and Gianluca Rinaldi
% Title: Why Does the Fed Move Markets so Much? 
% A Model of Monetary Policy and Time-Varying Risk Aversion (JFE, 2022)
%-------------------------------------------------------------------------%
clear all

tic
%% (Pflueger) Define the classes macro1, num1 and asset 
macro1  =   macro_dyn;
num1    =   num_set; 
asset   =   asset_p;

%% (Pfluger) Parameter settings that can speed up or slow down code 
% Define a dummy variable to run the risk neutral part (=1 for risk neutral)
asset.risk_neutral_run  =   0;

% Parameter settings for numerical solution method, e.g. grid density,
% number of simulations etc. 
num1 = num1.parameters;

% Fix random number generator if desired for exact replicability
rng('default');

% (Reichert) Define dummy to toggle lenghty simulation of irfs
simulate_irfs = 1;

%% (Reichert) Initialize Simulation results 

simulation_results          = struct();

%% (Pflueger) Input parameters for New Keynesian model 


macro1.rho_i            =   0.8;            % MP Persistence

% Preference parameters
macro1.theta0           =   0.9658;         % Peristence surplus consumption
macro1.theta1           =   -0.90;          % Backward looking habit
macro1.phi              =   0.9300;         % Consumption-output gap                                
macro1.gamma            =   2;              % Utility curvature
macro1.g                =   0.004725;       % Consumption growth      
macro1.rf               =   0.00235;        % Risk-free rate

% Phillips curve
macro1.kappa            =   0.0062/4;       % Slope of the Phillips curve
macro1.rho_pi           =   0.8;            % Backward-Looking Coefficient

% Leverage parameter
macro1.delta            =   0.6666; 

% Variance of quarterly MP shock in natural units. 
macro1.sigma_vec(3)     =   8.7767e-06 ;% 8.7767e-06 * 400    

% Std of FOMC date monetary policy surprise in annualized percent
asset.initialShockVec  = [0 0 0.0652 0]; 

% Compute implied Euler equation coefficients and fill in other parameters
% such as the industry portfolio betas
macro1 = macro1.update_params;

% Calibrate macroeconomic shock standard deviations
macro1.sigma_vec = [10^(-9) 10^(-9) macro1.sigma_vec(3) 10^(-9)]; % Base case from Pflueger2022

%% (Reichert) Set specification for gamma_x, gamma_y, for given calibration strategy

% Set dummy variable for simulation run with different levels of gamma_x
% and gamma_pi or the base parameters (=0)

% simulation_run = 0; Execute base case from Plueger2022
% simulation_run = 1; Execute for Strategy A or Sub-Strategies B1-B3, C1-C3
% simulation_run = 2; Execute for base case and Strategy A to plot IRFs
% simulation_run = 3; Execute for base case and Strategy B to plot IRFs
% simulation_run = 4; Execute for base case and Strategy C to plot IRFs
% simulation_run = 5; Execute for Strategy D: Grid of parameter values
simulation_run = 0;

%% (Reichert) Define list of possible parameters, while ensuring a unique solution
if simulation_run == 0
    gamma_x_list    = 0.5/4;
    gamma_pi_list   = 1.5;  
elseif simulation_run == 1 % use this to generate results for sub-strategies
    gamma_x_list    = 0;
    gamma_pi_list   = 1.001;
elseif simulation_run == 2
    gamma_x_list    = [0.5/4, 0];
    gamma_pi_list   = [1.5, 1.001 ];
elseif simulation_run == 3
    gamma_x_list    = [0.5/4, 1, 1.75, 2.50];
    gamma_pi_list   = [1.001, 1.5 ];
elseif simulation_run == 4
    gamma_x_list    = [0.5/4];
    gamma_pi_list   = [1.25, 1.5, 2.00, 2.75 ];
elseif simulation_run == 5
    gamma_x_list    = linspace(0, 2.5, 10);
    gamma_pi_list   = linspace(1.001, 2.5, 10); % starting point sufficiently larger than 1
end

%% (Reichert) Run Code for a set of parameters
for k = 1:length(gamma_pi_list)
    gamma_pi = gamma_pi_list(k);
    for j = 1:length(gamma_x_list)
        gamma_x = gamma_x_list(j);

        model_parameters = sprintf('gamma_x_%.1f_gamma_pi_%.1f', gamma_x, gamma_pi);
        model_parameters = strrep(model_parameters, '.', '_');
        
        % Monetary policy rule
        macro1.gamma_x          =   gamma_x;        
        macro1.gamma_pi         =   gamma_pi;       

    
        %% (Pflueger) Solve for macroeconomic dynamics of the form Y_t=PY_{t-1}+Qv_t
        
        % Solves for the macro dynamics of the model 
        macro1 = macro1.ModelPQ82(num1);
        
        %% (Pflueger) Prepare for asset price value function iteration
        
        % Compute the rotated state vector 
        macro1  = macro1.ScaledStateVector; 
        
        % Numerical settings for value function iterations, some of which depend on
        % rotated grid
        num1 = num1.update_num(macro1);
        
        %% (Pflueger) Solve and simulate risk-neutral asset prices.
        
        % Solve and simulate risk neutral asset prices. Only executed if
        % macro1.risk_neutral_run ==1
        %disp("Solve and simulate asset prices")
        tic
        asset = asset.risk_neutral_ap(macro1, num1); 
        toc
        %% (Pflueger) Solve and simulate full asset prices.
        
        % Implements the value function iteration and simulations
        disp('Computing prices')
        asset = asset.computeFn21(num1,macro1);   
        
        % (Reichert) Simulate path for macroeconomic dynamics and asset prices
        disp('Simulate moments')
        asset = asset.simulateMoments_wthout_FOMC(num1,macro1);

         %% (Pflueger) Structural Impulse responses. This computes the consumption moments in Table 3 (Pflueger)
        
        % Macroeconomic impulse responses do not need simulation. Set num.Nsim=2
        % and asset.testirf=1 for speed. 
        num1.Nsim = 2; 
        num1.Tirf = 150;
        asset.testirf=1;
        
        % Simulations to build impulse responses for asset prices
        asset =  asset.simulateStructuralIRF(macro1,num1);

        % Compare "Macro Irfs" with "Full-Scale Irfs"
        asset.Irf_saver = asset.Irftemp;

        % Lag trough consumption response
        lag_trough=find(asset.Irftemp.c==min(asset.Irftemp.c))-1;
        %% (Reichert) Save simulated asset moments in arrays to assess simultaniously (Reichert)

        %% Equity
        simulation_results.equity.eq_premium(k, j)      = asset.stocks.equityPremium;
        simulation_results.equity.vol(k, j)             = asset.stocks.vol;
        simulation_results.equity.sharpeRatio(k, j)     = asset.stocks.sharpeRatio;

        %% Nominal bonds
        simulation_results.nominal_bonds.term_premium(k,j)              = asset.nominalBonds.termPremium;
        simulation_results.nominal_bonds.vol(k, j)                      = asset.nominalBonds.vol;
        simulation_results.nominal_bonds.sharpeRatio(k, j)              = asset.nominalBonds.sharpeRatio;
        simulation_results.nominal_bonds.mean_log_yield_spread(k, j)    = asset.nominalBonds.meanLogYieldSpread;
        simulation_results.nominal_bonds.vol_log_yield_spread(k, j)     = asset.nominalBonds.volLogYieldSpread;
        simulation_results.nominal_bonds.coeffRegRetOnYS1y(k, j)        = asset.nominalBonds.coeffRegRetOnYS1y;
        
        %% Real Bonds
        simulation_results.real_bonds.term_premium(k,j)                 = asset.realBonds.termPremium;
        simulation_results.real_bonds.vol(k,j)                          = asset.realBonds.vol;
        simulation_results.real_bonds.sharpeRatio(k,j)                  = asset.realBonds.sharpeRatio;
        simulation_results.real_bonds.mean_log_yield_spread(k,j)        = asset.realBonds.meanLogYieldSpread;
        simulation_results.real_bonds.vol_log_yield_spread(k,j)         = asset.realBonds.volLogYieldSpread;
        
        %% Cross Asset
        simulation_results.CrossAsset.stock_bond_correlation(k,j)       = asset.crossAsset.corrNomStock;
        simulation_results.CrossAsset.beta_nom(k,j)                     = asset.crossAsset.betaNom;

        %% Macroeconomic Dynamics
        simulation_results.macro_dynamics.consGrowthVol = asset.macroDynamics.consGrowthVol;
        simulation_results.macro_dynamics.iChangeVol = asset.macroDynamics.iChangeVol;
        simulation_results.macro_dynamics.trough_output = min(100*asset.Irftemp.c)/max(asset.Irftemp.i);
        simulation_results.macro_dynamics.lag_trough = lag_trough;

        %% Macroeconomic Impulse Response Functions
        if k==1 && j == 1
            simulation_results.macro_irfs.x = [];
            simulation_results.macro_irfs.pi = [];
            simulation_results.macro_irfs.i = [];
            simulation_results.macro_irfs.pd = [];
            simulation_results.macro_irfs.y10nom = [];
        end
        simulation_results.macro_irfs.x = [simulation_results.macro_irfs.x, asset.Irftemp.x.'];
        simulation_results.macro_irfs.pi = [simulation_results.macro_irfs.pi, asset.Irftemp.pi.'];
        simulation_results.macro_irfs.i = [simulation_results.macro_irfs.i, asset.Irftemp.i.'];
        simulation_results.macro_irfs.pd = [simulation_results.macro_irfs.pd, asset.Irftemp.PD.'];
        simulation_results.macro_irfs.y10nom = [simulation_results.macro_irfs.y10nom, asset.Irftemp.y10nom.'];

        %% (Pflueger, adjusted by Reichert) Table 3: Quarterly moments
        if simulation_run == 0 || simulation_run == 1

            % Collect simulated model moments
            Table3      =   [asset.stocks.equityPremium, asset.stocks.vol, asset.stocks.sharpeRatio]';
            Table3      =   [Table3; 0;  [asset.nominalBonds.meanLogYieldSpread, asset.nominalBonds.vol, asset.nominalBonds.coeffRegRetOnYS1y]';...
                            0;asset.macroDynamics.consGrowthVol;asset.macroDynamics.iChangeVol];
            Table3      =   [Table3; 0; min(100*asset.Irftemp.c)/max(asset.Irftemp.i); lag_trough];
            
            % Data moments
            Table3Data  = [7.89; 
                          16.79; 
                          0.47; 
                          0; 
                          1.87;  
                          9.35;
                          2.69
                          0;
                          1.50;
                          1.35;
                          0;
                          -0.7;
                          4];
            
            % Create Table 3 to compare data and simulated model moments
            Table3 = [Table3, Table3Data];
            Table3 = array2table(round(Table3,2));
            Table3.Properties.RowNames = {'Equity Premium'...
            'Equity Vol'...
            'Equity SR'...
            '-'...
            'Yield Spread'...
            'Return Vol.'...
            '1-YR Excess REturns on Yield Spread'...
            '--'...
            'Std. Annual Cons. Growth'...
            'Std Annual Change Fed Funds Rate'...
            '---'...
            'Trough Effect Consumption'...
            'Lag Trough'
            };
            Table3.Properties.VariableNames = {'Model', 'Data'};
            Table3
        end
        
        %% (Reichert, Pflueger) Asset price impulse responses with simulations
        % (Pflueger) Set the parameters for simulations as the number of simulations and 
        % the length of IRFs
        num1.Nsim = 5000; 
        num1.Tirf = 150;
        asset.testirf=1;
        
        
        % (Reichert) Simulations to build impulse responses for asset prices
        if simulate_irfs == 1

            if k==1 && j==1
                asset_store_irfs = [];
            end 

            tic;
            disp("Simulate Impulse Response Functions")
            asset =  asset.simulateStructuralIRF(macro1,num1);
            toc;
            asset.Irf3=asset.Irftemp;
            
            asset_store_irfs.(model_parameters) = asset.Irf3;

            plot_StructuralIRF_STMP(asset)
        end
        

    end
end
toc

%% (Reichert) Plot asset moments given the monetary policy functions parameters 
if simulation_run == 2
    plot_strategy_a_irfs(simulation_results)
elseif simulation_run == 3
    plot_strategy_b_irfs(simulation_results)
elseif simulation_run == 4
    plot_strategy_c_irfs(simulation_results)
elseif simulation_run == 5
    plot_asset_moment_surfaces(simulation_results, gamma_x_list, gamma_pi_list)
end