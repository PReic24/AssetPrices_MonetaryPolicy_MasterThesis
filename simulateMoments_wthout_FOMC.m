%% (Reichert) Function simulating asset moments, without including FOMC specific calculations, 
%% Otherwise identical to simulateMoments in asset_p.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function simulateMoments_wthout_FOMC(asset_pI, num_set, macro_dyn)
            
            asset_pI.risk_neutral_run   = asset_pI.risk_neutral_run_aux;
            Nsim    = num_set.Nsim;
            T       = num_set.T;
            N       = num_set.N;
            sizexm  = num_set.sizexm;
            shat    = num_set.shat;
            xmgrid  = num_set.xmgrid;
            N2      = num_set.N2;
            burn    = num_set.burn;
            z       = num_set.Z;
            Nbonds  = num_set.Nbonds;
            Nbins   = num_set.Nbins;
            h       = num_set.h;
            gamma     = macro_dyn.gamma;
            delta     = macro_dyn.delta;
            if isempty(macro_dyn.delta_betaport)
                % If the leverage of beta-sorted portfolios is uninitialized,
                % then set it to be the leverage of the market portfolio.
                delta_betaport = macro_dyn.delta;
            else
                delta_betaport = macro_dyn.delta_betaport;
            end
            theta0    = macro_dyn.theta0;
            P         = macro_dyn.P;
            Q         = macro_dyn.Q;
            A         = macro_dyn.A;
            Ainv      = macro_dyn.Ainv;
            Sigmau    = macro_dyn.Sigmau;
            Ptilde    = macro_dyn.Ptilde;
            g         = macro_dyn.g;
            phi       = macro_dyn.phi;
            sigmac    = macro_dyn.sigmac;
            rf        = macro_dyn.rf;
            Sbar      = macro_dyn.Sbar;
            smax      = macro_dyn.smax;
            QM        =[1,0,0]*macro_dyn.Q;
            sigma_vec = macro_dyn.sigma_vec;
            initialShockVec = asset_pI.initialShockVec;
            
            % Pre-allocate constants used for dynamics of shat
            const2 =  [1,0,0]*( P- phi*eye(3))* Ainv;
            const3 = ([0,0,1]-[0,1,0]* P)* Ainv;
            
            % Copy asset prices
            asset = asset_p(asset_pI.G,asset_pI.Bn,asset_pI.Bnom,macro_dyn.ImpliedParams);
            asset.G_rn              = asset_pI.G_rn;
            asset.Bn_rn             = asset_pI.Bn_rn;
            asset.Bnom_rn           = asset_pI.Bnom_rn;
            asset.rho               = asset_pI.rho;
            asset.FOMC_run          = asset_pI.FOMC_run;
            asset.risk_neutral_run  = asset_pI.risk_neutral_run;
            asset.simulated_rn      = asset_pI.simulated_rn;
                       
            % Points and weights for numerical integration over MP shock
            num_set            = num_set.generateprobs;
            xGL                = num_set.xGL2;
            
            % Define matrices used often
            Z14                = zeros(1,4);
            Z3T                = zeros(3,T);
            ZT1                = zeros(T,1);                     
            
            % Reshape price-consumption ratio
            G5dim              = log(reshape(asset_pI.G, N, N, N,N2,sizexm));
            
            % Reshape nominal bond prices
            % Horizon for RHS yields in real information effect regressions
          
            Pnomplus5dim       = log(reshape(asset_pI.Bnom(end,:,:,:), N ,N ,N, N2, sizexm));
            Pnomminus5dim      = log(reshape(asset_pI.Bnom(end-1,:,:,:), N ,N ,N, N2, sizexm));
            Pnomhq5dim         = log(reshape(asset_pI.Bnom(h-1,:,:,:), N ,N ,N, N2, sizexm));
            Pnom5y5dim         = log(reshape(asset_pI.Bnom(19,:,:,:), N ,N ,N, N2, sizexm));
            
            % Reshape real bond prices
            Pplus5dim          = log(reshape(asset_pI.Bn(end,:,:,:), N, N, N, N2, sizexm));
            Pminus5dim         = log(reshape(asset_pI.Bn(end-1,:,:,:), N, N, N, N2, sizexm));
            P5y5dim            = log(reshape(asset_pI.Bn(19,:,:,:), N ,N ,N, N2, sizexm));
            
           
            % Do the same for risk neutral asset prices. If no rn then just give the same
            
            if asset.risk_neutral_run == 1
                
                G5dim_rn            = log(reshape(asset_pI.G_rn, N, N, N,N2,sizexm));
                
                Pnomplus5dim_rn     = log(reshape(asset_pI.Bnom_rn(end,:,:,:), N ,N ,N, N2, sizexm));
                Pnomminus5dim_rn    = log(reshape(asset_pI.Bnom_rn(end-1,:,:,:), N ,N ,N, N2, sizexm));
                Pnomhq5dim_rn       = log(reshape(asset_pI.Bnom_rn(h-1,:,:,:), N ,N ,N, N2, sizexm));
                Pnom5y5dim_rn       = log(reshape(asset_pI.Bnom_rn(19,:,:,:), N ,N ,N, N2, sizexm));
                
                Pplus5dim_rn        = log(reshape(asset_pI.Bn_rn(end,:,:,:), N, N, N, N2, sizexm));
                Pminus5dim_rn       = log(reshape(asset_pI.Bn_rn(end-1,:,:,:), N, N, N, N2, sizexm));
                P5y5dim_rn          = log(reshape(asset_pI.Bn_rn(19,:,:,:), N ,N ,N, N2, sizexm));
            else
                
                G5dim_rn             = log(reshape(asset_pI.G, N, N, N,N2,sizexm));
                
                Pnomplus5dim_rn      = log(reshape(asset_pI.Bnom(end,:,:,:), N ,N ,N, N2, sizexm));
                Pnomminus5dim_rn     = log(reshape(asset_pI.Bnom(end-1,:,:,:), N ,N ,N, N2, sizexm));
                Pnomhq5dim_rn        = log(reshape(asset_pI.Bnom(h-1,:,:,:), N ,N ,N, N2, sizexm));
                Pnom5y5dim_rn        = log(reshape(asset_pI.Bnom(19,:,:,:), N ,N ,N, N2, sizexm));
                
                Pplus5dim_rn         = log(reshape(asset_pI.Bn(end,:,:,:), N, N, N, N2, sizexm));
                Pminus5dim_rn        = log(reshape(asset_pI.Bn(end-1,:,:,:), N, N, N, N2, sizexm));
                P5y5dim_rn           = log(reshape(asset_pI.Bn(19,:,:,:), N ,N ,N, N2, sizexm));
            end
            
            
            % Reshape grid for n-dimensional interpolation
            z1shape            = reshape(z(:,1), N, N, N);
            z2shape            = reshape(z(:,2), N, N, N);
            z3shape            = reshape(z(:,3), N, N, N);
            z1grid             = reshape(z1shape(:,1,1), N, 1);
            z2grid             = reshape(z2shape(1,:,1), N, 1);
            z3grid             = reshape(z3shape(1,1,:), N, 1);
            
            % Upper and lower bounds of grids
            zgrid_             = [z1grid, z2grid, z3grid];
            zupper             = max(zgrid_)';
            zlower             = min(zgrid_)';
            supper             = max(shat);
            slower             = min(shat);
            xminusupper        = max(xmgrid);
            xminuslower        = min(xmgrid);
            
            % Start Nsim independent simulations
            rng(0) % To avoid simulation noise
            % j counts the independent simulations
            j=1;
            while j<=1    
                %% Draw shocks - we split all shocks into the component before and after the MP date
                % Generate T draws of u_t and compute \epsilon_t from that
                u            = mvnrnd(Z14, diag(sigma_vec), T)';
                % Information released on MP date
                eps      = A*Q*u;
                uast     = u(4,:);
                                
                % Initialize simulated time series, time series with suffix '_pre' denote prices just prior to MP release
                ztildesim      = Z3T;
                shatsim        = ZT1;
                shatsim_pre    = ZT1;
                Yhat           = Z3T;
                Yhat_pre       = Z3T;
                PD             = ZT1;
                PD_pre         = ZT1;
                PD_rn_pre      = ZT1;
                PD_rn          = ZT1;
                Pnomplus       = ZT1;
                Pnomplus_pre   = ZT1;
                Pnomplus_rn_pre= ZT1;
                Pnomplus_rn    = ZT1;
                Pnomminus      = ZT1;
                Pnomminus_rn   = ZT1;
                Pplus          = ZT1;
                Pplus_pre      = ZT1;
                Pplus_rn_pre   = ZT1;
                Pplus_rn       = ZT1;
                Pminus         = ZT1;
                Pminus_rn      = ZT1;
                csim           = ZT1;
                csim_pre       = ZT1;
                piast          = ZT1;
                piast_pre      = ZT1;
                approxErrI     = ZT1;
                Pnomhq         = ZT1;
                Pnomhq_rn      = ZT1;
                Pnomhq_pre     = ZT1;
                Pnomhq_rn_pre  = ZT1;
                P5y            = ZT1;
                P5y_rn         = ZT1;
                P5y_pre        = ZT1;
                P5y_rn_pre     = ZT1;
                Pnom5y         = ZT1;
                Pnom5y_rn      = ZT1;
                Pnom5y_pre     = ZT1;
                Pnom5y_rn_pre  = ZT1;
                
                % Update state vector
                for t=3:T
                    % Dynamics for \tilde Z
                    ztildesim(:,t) = Ptilde*ztildesim(:,t-1)+eps(:,t);
                    
                    % Dynamics for \hat Y
                    Yhat(:,t)      = P*Yhat(:,t-1)+Ainv*eps(:,t);
                    % Dynamics for surplus consumption ratio relative to steady-state
                    shatsim(t)     = ...
                        theta0*shatsim(t-1)+((1/gamma)*const3-const2)*ztildesim(:,t-1)+...
                        senshat(shatsim(t-1), Sbar)*sigmac*eps(1,t);                    
                    % Truncate state variables at upper and lower end of grid, so we can use standard interpolation
                    zinterp        = max(zlower, min(ztildesim(:,t), zupper));
                    sinterp        = max(slower, min(shatsim(t), supper));
                    xminusinterp   = max(xminuslower, min(const3*ztildesim(:,t-1),xminusupper)); 
                    csim(t)        = g+csim(t-1)+(Yhat(1,t)-phi*Yhat(1,t-1));
                    piast(t)       = piast(t-1)+uast(t);
                    
                    % Scaling factors for nominal bonds
                    ePplus    = exp(Nbonds*piast(t));
                    ePminus   = exp((Nbonds-1)*piast(t));
                    eP5y      = exp(20*piast(t));
                    ePhq      = exp(h*piast(t));
                    
                    % Interpolate to obtain price-consumption ratio and real and nominal bond prices with risk premia
                    qNDInterp = eqinterpnND2(z1grid, z2grid, z3grid, shat, xmgrid, ...
                        {G5dim, Pnomplus5dim, Pnomminus5dim, Pplus5dim, Pminus5dim, Pnomhq5dim, P5y5dim, Pnom5y5dim}, ...
                        zinterp(1), zinterp(2), zinterp(3), ...
                        sinterp, xminusinterp, 0);
                    
                    % Scale price-consumption ratio for consumption claim
                    PD(t)          = qNDInterp(1)/4;
                    
                    % n-quarter bond prices
                    Pnomplus(t)    = qNDInterp(2)/ePplus;
                    Pplus(t)       = qNDInterp(4);
                    
                    % n-1-quarter bond prices
                    Pnomminus(t)   = qNDInterp(3)/ePminus;
                    Pminus(t)      = qNDInterp(5);
                    
                    % h quarters bond prices
                    Pnomhq(t) = qNDInterp(6)/ePhq;
                    
                    % 5 year bond prices
                    P5y(t) = qNDInterp(7);
                    Pnom5y(t) = qNDInterp(8)/eP5y;
                    
                    % Interpolate to obtain risk neutral price-consumption ratio and real and nominal bond prices
                    qNDInterp_rn = eqinterpnND2(z1grid, z2grid, z3grid, shat, xmgrid, ...
                        {G5dim_rn, Pnomplus5dim_rn, Pnomminus5dim_rn, Pplus5dim_rn, Pminus5dim_rn, Pnomhq5dim_rn, P5y5dim_rn, Pnom5y5dim_rn}, ...
                        zinterp(1), zinterp(2), zinterp(3), ...
                        sinterp, xminusinterp, 0);
                    
                    %scale risk-neutral price-consumption ratio for consumption claim
                    PD_rn(t)          = qNDInterp_rn(1)/4;
                    
                    %risk-neutral n-quarter bond prices
                    Pnomplus_rn(t)    = qNDInterp_rn(2)/ePplus;
                    Pplus_rn(t)       = qNDInterp_rn(4);
                    
                    %risk-neutral n-1-quarter bond prices
                    Pminus_rn(t)      = qNDInterp_rn(5);
                    Pnomminus_rn(t)   = qNDInterp_rn(3)/ePminus;
                    
                    % h quarters risk neutral bond prices
                    Pnomhq_rn(t) = qNDInterp_rn(6)/ePhq;
                    
                    % 5 year risk neutral bond prices
                    P5y_rn(t) = qNDInterp_rn(7);
                    Pnom5y_rn(t) = qNDInterp_rn(8)/eP5y;
                    
                    % approximation error for 1-period nominal rate
                    approxErrI(t) = .5*([0,1,0]*Q(:,2:4) + [0,0,1])*Sigmau(2:4,2:4)*([0,1,0]*Q(:,2:4) + [0,0,1])'+ ...
                        gamma*(senshat(shatsim(t), Sbar) + 1)*QM(2:4)*Sigmau(2:4,2:4)*([0,1,0]*Q(:,2:4) + [0,0,1])';
                    
                end
                % Nbonds/4-year nominal and real bond yields in annualized percent (in logs!)
                y10nom              = -400*log(Pnomplus)/Nbonds;
                y10real             = -400*log(Pplus)/Nbonds;
                y10nom_rn           = -400*log(Pnomplus_rn)/Nbonds;
                y10real_rn          = -400*log(Pplus_rn)/Nbonds;
                % h quarter nominal bonds yields
                yhqnom              = -400*log(Pnomhq)/h;
                yhqnom_rn           = -400*log(Pnomhq_rn)/h;
                % 5 year bond yields
                y5nom               = -400*log(Pnom5y)/20;
                y5nom_rn            = -400*log(Pnom5y_rn)/20;
                y5real              = -400*log(P5y)/20;
                y5real_rn           = -400*log(P5y_rn)/20;
                
                % Real short rate (i minus expected inflation)
                rfr                = 400*([0,0,1]*Yhat-[0,1,0]*P*Yhat+rf);
                
                % Levered dividends: Compute D^{delta}_{t+1}
                dDelta = PD(2:end).*exp(csim(2:end)) + exp(csim(2:end)) - (1-delta).*PD(1:end-1).*exp(csim(1:end-1)).*exp(rfr(1:end-1)./400)' ...
                    - delta.*PD(2:end).*exp(csim(2:end));
                
                % Take the 64-quarter moving average of D^{delta}_{t+1}
                dDelta_bar = conv(dDelta,ones(1,64),'valid')/64;
                
                if burn>200
                    dDelta_bar = 0.5*(dDelta_bar(1:end-64)+dDelta_bar(65:end));
                end

                % Drop burn period observations and other observations to ensure that returns and prices have
                % the same length
                PD                 = PD(burn-1+2:end);
                PD_rn              = PD_rn(burn-1+2:end);
                Yhat               = Yhat(:,burn-1+2:end);
                shatsim            = shatsim(burn-1+2:end);
                Pnomplus           = Pnomplus(burn-1+2:end);
                Pnomminus          = Pnomminus(burn-1+2:end);
                Pplus              = Pplus(burn-1+2:end);
                Pminus             = Pminus(burn-1+2:end);
                Pnomplus_rn        = Pnomplus_rn(burn-1+2:end);
                Pnomminus_rn       = Pnomminus_rn(burn-1+2:end);
                Pplus_rn           = Pplus_rn(burn-1+2:end);
                Pminus_rn          = Pminus_rn(burn-1+2:end);
                csim               = csim(burn-1+2:end);
                piast              = piast(burn-1+2:end)';
                y10nom             = y10nom(burn-1+2:end);
                y10real            = y10real(burn-1+2:end);
                rfr                = rfr(burn-1+2:end);
                eps                = eps(:,burn-1+3:end);
                approxErrI         = approxErrI(burn-1+2:end); 
                burn1              = burn+2; 

                if burn>200
                    dDelta_bar         = dDelta_bar((burn-129)+2:end); 
                else
                    dDelta_bar         = dDelta_bar((burn-65)+2:end); 
                end
                
                % Price-dividend ratio of levered equity at time t divides by smoothed dividends
                PDlev              = delta.*PD.*exp(csim)./dDelta_bar;
                
                % 1-quarter nominal yield
                rfr_nom            = 400*([0,0,1]*Yhat+piast+rf);
                
                % Level return on consumption claim
                Ret                = exp((csim(2:end)-csim(1:end-1))).*(1+4*PD(2:end))./(4*PD(1:end-1));
                Ret_rn             = exp((csim(2:end)-csim(1:end-1))).*(1+4*PD_rn(2:end))./(4*PD_rn(1:end-1));
                
                % Log return on levered equity
                reteq              = log((1/delta)*Ret - ((1-delta)/delta)*exp(rfr(1:end-1)'/400));
                reteq_rn           = log((1/delta)*Ret_rn - ((1-delta)/delta)*exp(rfr(1:end-1)'/400));
                delta_temp = repmat(delta_betaport', length(Ret), 1);
                Ret_temp = repmat(Ret, 1, length(delta_betaport));
                rfr_temp = repmat(rfr(1:end-1)'/400, 1, length(delta_betaport));
                ret_betaport = log((1./delta_temp).*Ret_temp - ((1-delta_temp)./delta_temp).*exp(rfr_temp)); 
                
                % In percentage units in excess of the riskfree rate
                reteq              = 100*reteq - rfr(1:end-1)'/4;
                reteq_rn           = 100*reteq_rn - rfr(1:end-1)'/4;
                ret_betaport = 100*ret_betaport - rfr(1:end-1)'/4; %XXXFL
                
                % Log excess nominal bond returns in percent but not annualized
                retnom             = 100*log(Pnomminus(2:end))-100*log(Pnomplus(1:end-1))-rfr_nom(1:end-1)'/4;
                retnom_rn          = 100*log(Pnomminus_rn(2:end))-100*log(Pnomplus_rn(1:end-1))-rfr_nom(1:end-1)'/4;
                
                % Log excess real bond returns in percent but not annualized
                retreal            = 100*log(Pminus(2:end))-100*log(Pplus(1:end-1))-rfr(1:end-1)'/4;
                retreal_rn         = 100*log(Pminus_rn(2:end))-100*log(Pplus_rn(1:end-1))-rfr(1:end-1)'/4;
                
                % Breakeven: log returns on nominal in excess of log returns on real bonds
                breakeven          = retnom - retreal;

                % Nominal and real log yield spreads
                spreadNom          	= y10nom'-rfr_nom;
                spreadReal          = y10real'-rfr;

                % 1-year log equity excess returns in natural units
                ret1yr             = conv(reteq,ones(1,4),'valid')/100;
                
                % Levered price-dividend ratio
                pdlev = log(PDlev);
                
                % Compute cash-flow news and real rate news vectors to compute equity real rate news analytically according to Campbell and Ammer
                rho                = asset_pI.rho;
                rhoIPinv           = inv(eye(3)-rho*P);
                Gammaeq_rr         = -rho*([0,0,1]-[0,1,0]*P)*rhoIPinv*Ainv;
                
                % Cash-flow news of stock and bond returns
                reteq_cf           = reteq_rn-(100*Gammaeq_rr*eps)';
                retnom_cf          = retnom_rn-retreal_rn;
                
                % Risk premium excess returns of stock and bond returns
                reteq_rp           = reteq-reteq_rn;
                retnom_rp          = retnom-retnom_rn;
                                            
                                               
                %% Stock moments
                
                % Equity risk premium in annualized units
                EqPremium       = 4*(mean(reteq) + .5*std(reteq)^2/100);
                
                % Std. log equity excess returns
                Stdeq           = std(reteq)*2;
                
                % Exp(mean(log pd))
                mean_pdlev      = exp(mean(pdlev));
                
                % Standard deviation of log dp
                std_dp          = std(pdlev);
                
                % Autocorrelation of dp
                dp_corr            = corrcoef(pdlev(2:end), pdlev(1:end-1));
                rho_dp          = dp_corr(1,2);
                
                % Predictability with pd 1 quarter, 1 year and 5 year regressions of returrn on price-dividend ratios
                [re1_coef,~,~,~,R2_re1] = regress(ret1yr, [ones(size(pdlev(1:end-4))), pdlev(1:end-4)]);
                re1                  = re1_coef(2);
                re1_r2               = R2_re1(1);
                
                %% Nominal bond moments
                % Nominal term premium
                BondPremium     = 4*(mean(retnom) + .5*std(retnom)^2/100);
                % Std of bond returns
                Stdnom          = std(retnom)*2;
                
                % Mean and std of log yield spread
                TermSlope       = mean(spreadNom);
                TermSlopeStd    = std(spreadNom);
                
                % Autocorrelation of log yield spread
                AR_slope5temp      = corrcoef(spreadNom(2:end), spreadNom(1:end-1));
                AR_slope5       = AR_slope5temp(1,2);
                
                % Regress returns onto lagged log yield spread
                ret1yrNom          = retnom(4:end)+retnom(3:end-1)+retnom(2:end-2)+retnom(1:end-3);
                ret1yrNom          = ret1yrNom/100;
                % Multiply returns by 100 to match units in empirical exercise
                [ys1_coef,~,~,~,R2_ys1] = regress(100*ret1yrNom, [ones(size(spreadNom(1:end-4)')), spreadNom(1:end-4)']);
                ys1             = ys1_coef(2);
                ys1_r2          = R2_ys1(1);
                
                %bond returns are predicted by PD ratio
                [ypd1_coef,~,~,~,R2_ypd1] = regress(100*ret1yrNom, [ones(size(spreadNom(1:end-4)')), pdlev(1:end-4)]);
                ypd1             = ypd1_coef(2);
                ypd1_r2          = R2_ypd1(1);
                
                [ypd2_coef,~,~,~,R2_ypd2] = regress(100*(ret1yrNom(1:end-4)+ret1yrNom(5:end)), [ones(size(spreadNom(1:end-8)')), pdlev(1:end-8)]);
                ypd2             = ypd2_coef(2);
                ypd2_r2          = R2_ypd2(1);
                
                %% Cross-asset
                % Bond-stock return correlations
                bondstock_corr_temp     = corrcoef(retnom, reteq);
                tipsstock_corr_temp     = corrcoef(retreal, reteq);
                correlations       = [bondstock_corr_temp(1,2), tipsstock_corr_temp(1,2)];
                
                % Nominal bond beta
                beta_temp           = regress(retnom, [ones(T-burn1+1,1), reteq]);
                beta_nom         = beta_temp(2);
                
                % 1-year excess stock return on output gap
                [coeffStockGap_temp,~,~,~,R2StockGap_temp]         = regress(ret1yr, [ones(size(Yhat(1,1:end-4)')), Yhat(1,1:end-4)']);
                coeffStockGap   = coeffStockGap_temp(2);
                R2StockGap      = R2StockGap_temp(1);
                
                % 1-year excess bond return on output gap
                [coeffBondGap_temp,~,~,~,R2BondGap_temp]         = regress(ret1yrNom, [ones(size(Yhat(1,1:end-4)')), Yhat(1,1:end-4)']);
                coeffBondGap    = coeffBondGap_temp(2);
                R2BondGap       = R2BondGap_temp(1);
                
                % Betas of beta-sorted portfolios.
                beta_temp = cell2mat(cellfun(@(y) regress(y, [ones(length(reteq),1), reteq]), ...
                    num2cell(ret_betaport, 1), ...
                    'UniformOutput', false));
                beta_port = beta_temp(2, :)';
                
                %% Real bonds moments
                % Term premium
                RealBondPremium = 4*(mean(retreal) + .5*std(retreal)^2/100);
                % Std returns
                Stdreal         = std(retreal)*2;
                
                % Mean and std. log yield spread
                TermSlopeReal   = mean(spreadReal);
                TermSlopeRealStd= std(spreadReal);
                % Real bond beta
                beta_temp          = regress(retreal, [ones(T-burn1+1,1), reteq]);
                beta_real       = beta_temp(2);
                
                %% Breakeven moments
                BreakevenTermPremium      =  4*(mean(breakeven) + .5*std(breakeven)^2/100);
                StdBreakeven              =  std(breakeven)*2;
                SharpeRatioBreakeven      =  BreakevenTermPremium(j)/StdBreakeven(j);
                
                beta_temp                    = regress(breakeven, [ones(T-burn1+1,1), reteq]);
                beta_breakeven            = beta_temp(2);
                breakevenSim(1:(T-burn1+2))  = y10nom - y10real;
                %% Macro dynamics
                % Std and AR(1) of changes in nominal 1-quarter yield
                rfrNomStd       = std(rfr_nom(5:end)-rfr_nom(1:end-4));
                
                % Std inflation changes
                piChanges            = 4*(100*Yhat(2,1:end-4)+100*piast(1:end-4) - (100*Yhat(2,5:end)+100*piast(5:end)));
                piChangeVol       = std(piChanges);
                         
                % Std Log consumption growth
                consGrowthVol     = 100*std(csim(5:end)-csim(1:end-4));
                
                % Std Output gap
                xVol              = 100*std(Yhat(1,5:end)-Yhat(1,1:end-4));
                
                % Slope coefficient consumption growth output gap growth
                slopetemp            = regress((csim(5:end)-csim(1:end-4)), [ones(size(Yhat(1,5:end))); Yhat(1,5:end)-Yhat(1,1:end-4)]');
                slope_cx             = slopetemp(2,1);

                % Correlation(Delta rfr, Delta consumption growth)
                corrtemp                 = corrcoef(csim(5:end)'-csim(1:end-4)', rfr_nom(5:end)-rfr_nom(1:end-4));
                corr_c_rfr            = corrtemp(1,2);
                
                % Fraction s_t>s^max
                fracmax         = sum(shatsim+log(Sbar)>smax)/T;
                % Std of approximation error for 1-period nominal rate in annualized percent
                stdApproxErrI     = std(approxErrI)*400;
                
                % Skip simulation run if PDlev turns negative
                if min(PDlev)<0
                    j=j+1;
                    continue
                end
                % Increase loop counter
                j=j+1;
            end
            %% Save output
            asset.additional_moments     = [mean(fracmax), stdApproxErrI];
            
            %% Equities
            stocks.equityPremium     = EqPremium;
            stocks.vol               = mean(Stdeq);
            stocks.sharpeRatio       = EqPremium/mean(Stdeq);
            stocks.meanPDlev         = mean_pdlev;
            stocks.stdDP             = std_dp;
            stocks.rhoDP             = rho_dp;
            stocks.coeffRegRetOnPD1y = re1;
            stocks.R2RegRetOnPD1y    = re1_r2;
            
            asset.stocks             = stocks;
            
            %% Nominal Bonds
            nominalBonds.termPremium               = BondPremium;
            nominalBonds.vol                       = mean(Stdnom);
            nominalBonds.sharpeRatio               = BondPremium/mean(Stdnom);
            nominalBonds.meanLogYieldSpread        = TermSlope;
            nominalBonds.volLogYieldSpread         = TermSlopeStd;
            nominalBonds.persistenceLogYieldSpread = AR_slope5;
            nominalBonds.betaNom                   = beta_nom;
            nominalBonds.coeffRegRetOnYS1y         = ys1;
            nominalBonds.R2RegRetOnYS1y            = ys1_r2;
            nominalBonds.coeffRegRetOnPD1y         = ypd1;
            nominalBonds.R2RegRetOnPD1y            = ypd1_r2;
            nominalBonds.coeffRegRetOnPD2y         = ypd2;
            nominalBonds.R2RegRetOnPD2y            = ypd2_r2;
            
            asset.nominalBonds                = nominalBonds;
            
            %% Real Bonds
            realBonds.termPremium          = RealBondPremium;
            realBonds.vol                  = mean(Stdreal);
            realBonds.sharpeRatio          = RealBondPremium/mean(Stdreal);
            realBonds.meanLogYieldSpread   = TermSlopeReal;
            realBonds.volLogYieldSpread    = TermSlopeRealStd;
            realBonds.betaRealStock        = beta_real;
            realBonds.corrRealStock        = correlations(2);
            
            asset.realBonds = realBonds;
            
            %% Brekaven moments
            breakevens.termPremium          = BreakevenTermPremium;
            breakevens.vol                  = StdBreakeven;
            breakevens.sharpeRatio          = SharpeRatioBreakeven;
            breakevens.stockBeta            = beta_breakeven;
            breakevens.simulation           = breakevenSim;
            asset.breakevens = breakevens;
            
            %% Cross-Asset moments
            crossAsset.corrNomStock     = correlations(1);
            crossAsset.betaNom          = beta_nom;
            crossAsset.coeffStockGap    = coeffStockGap;
            crossAsset.R2StockGap       = R2StockGap;
            crossAsset.coeffBondGap     = coeffBondGap;
            crossAsset.R2BondGap        = R2BondGap;
            crossAsset.beta_port        = beta_port; %XXXFL
            
            asset.crossAsset = crossAsset;
            
            %% Macro Dynamics moments
            macroDynamics.iChangeVol     = rfrNomStd;
            macroDynamics.corr_c_rfr     = corr_c_rfr;    
            macroDynamics.piChangeVol    = mean(piChangeVol);
            macroDynamics.consGrowthVol  = mean(consGrowthVol);
            macroDynamics.xVol           = mean(xVol);
            macroDynamics.slope_cx       = mean(slope_cx);
            
            asset.macroDynamics = macroDynamics;