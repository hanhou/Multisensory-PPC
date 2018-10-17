function [fitted_w, cost] = LinearFit(fit_optimal_trace_interp,fit_basis_trace_per_cell,alpha)
% Do fitting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if_plot = 0;
swap_visual_vest = 0; % Two fellows of Alex suggested us try flip the vest and vis labels HH20180908
increase_step_trick = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

w0 = increase_step_trick + (1+randn(1,size(fit_basis_trace_per_cell,1)))/size(fit_basis_trace_per_cell,1); % Start from straight average
proj_PSTH = reshape(w0 * reshape(fit_basis_trace_per_cell,size(fit_basis_trace_per_cell,1),[]),[],3);

if if_plot
    opt = optimset('MaxFunEvals',50000,'MaxIter',300,'Display','off','PlotFcns',@f7p2_fitting_deltaPSTH_plot_func);
    set(figure(1458),'name','Optimization PlotFcns'); clf;
else
    opt = optimset('MaxFunEvals',50000,'MaxIter',300,'Display','off','PlotFcns',[]);
end

if swap_visual_vest == 0 % Normal labeling
    
    [fitted_w, cost, exitflag, output] = fmincon(@(w)CostFunc(fit_optimal_trace_interp,...
        fit_basis_trace_per_cell,w,alpha,increase_step_trick),...
        w0,[],[],[],[],0*w0,[],[],opt);
    
else % Total flip, and only fit vest and vis
    ys = fit_optimal_trace_interp;
    ys(:,3) = nan; % Don't fit combined
    
    xs = fit_basis_trace_per_cell(:,:,[2,1,3]);
    xs(:,:,3) = nan;
    
    [fitted_w, cost, exitflag, output] = fmincon(@(w)CostFunc(ys,xs,w,alpha,increase_step_trick),...
        w0,[],[],[],[],0*w0,[],[],opt);
end

%     scan_alpha_costs{aa} = fitting_dynamics{1};
%     scan_alpha_weight{aa} = fitted_w;


    function stop = f7p2_fitting_deltaPSTH_plot_func(weight,optimValue,~)
        %global fitting_dynamics increase_step_trick;
        
        weight = weight-increase_step_trick;
        
        persistent costs;
        if optimValue.iteration == 0
            costs = [];
        end
        set(findobj(gcf,'type','axes'),'visible','off');
        h = tight_subplot(2,2,[0.1 0.1],0.1,0.1);
        
        % Update cost
        proj_test_PSTH = reshape((weight) * reshape(test_basis_trace_per_cell,size(test_basis_trace_per_cell,1),[]),[],3);
        test_cost = mean((proj_test_PSTH(:) - test_optimal_trace_interp(:)).^2); % Mean Squared error
        proj_valid_PSTH = reshape((weight) * reshape(valid_basis_trace_per_cell,size(valid_basis_trace_per_cell,1),[]),[],3);
        valid_cost = mean((proj_valid_PSTH(:) - valid_optimal_trace_interp(:)).^2); % Mean Squared error
        costs = [costs; optimValue.fval-alpha*norm(weight) test_cost valid_cost];
        
        if mod(optimValue.iteration,1) == 0
            % Fitted delta PSTH
            %             plot(h(1),real_ts,model_PSTH_interp,'linew',2); hold(h(1),'on');
            %             plot(h(1),real_ts,reshape((weight) * reshape(real_deltaPSTHs_per_cell,size(real_deltaPSTHs_per_cell,1),[]),[],3),'--','linew',2);
            plot(h(1),data_to_fit_time,interp1(model_ts,model_PSTH_optimal_M1(:,:),data_to_fit_time),'linew',2); hold(h(1),'on');
            plot(h(1),data_to_fit_time,reshape((weight) * reshape(data_to_fit_PSTH,size(data_to_fit_PSTH,1),[]),[],3),'--','linew',2);
            xlim(h(1),[-400 2400]);
            ylim(h(1),[-5 20]);
            ylims = ylim((h(1))); indicator_pos = min(ylims)*ones(1,length(common_ts));
            plot(h(1),[min(model_ts) min(model_ts)],[min(ylims) max(ylims)],'k-');
            plot(h(1),[max(model_ts) max(model_ts)],[min(ylims) max(ylims)],'k-');
            plot(h(1),common_ts(fit_time_range),indicator_pos(fit_time_range),'ok','markerfacecol','k','markersize',5);
            plot(h(1),common_ts(test_time_range),indicator_pos(test_time_range),'or','markerfacecol',colors(2,:),'markersize',5);
            plot(h(1),common_ts(valid_time_range),indicator_pos(valid_time_range),'ob','markerfacecol',colors(1,:),'markersize',5);
            % plot(h(1),[min(real_ts) max(real_ts)],min(ylims)*ones(1,2),'k-','linew',5);
            title(h(1),sprintf('Solid: model; Dashed: %s; Swap = %g',data_to_fit{use_data_to_fit,3}, swap_visual_vest));
            ylabel(h(1),'delta firing');
            
            
            % Fitting error
            plot(h(3), 1:size(costs,1), log(costs(:,1)),'ok-');  % Training error
            hold(h(3),'on');
            plot(h(3), 1:size(costs,1), log(costs(:,2)),'or-');  % Testing error
            plot(h(3), 1:size(costs,1), log(costs(:,3)),'ob-');  % Testing error
            
            legend(h(3),{'Training','Testing','Validating'});
            title(h(3),sprintf('\\alpha = %g, step = %g',alpha,optimValue.iteration));
            ylabel(h(3),'log(mean of squared error)');
            xlabel(h(3),'# Iteration');
            
            % Weight
            plot(h(4),sort(weight),'o'); hold(h(4),'on');
            xlabel(h(4),'# cell');
            ylabel(h(4),'Fitting weight');
            
        end
        stop = false;
        
        % Save data
        fitting_dynamics = {costs};
        
    end

end