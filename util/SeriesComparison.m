function result = SeriesComparison(ys,ts,varargin)
% Compare time series result = SeriesComparison(ys,ts,varargin) 
% If ys is matrix, the third dimension represents the category HH20141224
%
% ------ Input parameters -------
% paras = inputParser;
% 
% addOptional(paras,'OverrideError',[]); % Error overrided by inputs. HH20160908
% addOptional(paras,'OverridePs',[]); % P values overrided by inputs. HH20160908
% 
% addOptional(paras,'figN',999);
% addOptional(paras,'axes',[]);
% 
% addOptional(paras,'ErrorBar',2);  % Sum of (1: Normal Errorbar; 2: Shaded Errorbar; 4: Significance marker)
% addOptional(paras,'CompareIndex',[]);  % Pair of statistic test for drawing significance marker: {[1 3 5; 2 4 6]}. HH20160427
% addOptional(paras,'CompareColor',{});  % Colors for significance marker.
% addOptional(paras,'PCritical',0.01);  % Colors for significance marker.
% 
% addOptional(paras,'Colors',{'b','r','g'}); % (originally for comparing PSTH of LIP neurons)
% addOptional(paras,'LineStyles',{'-'});
% addOptional(paras,'Transparent',1);
% addOptional(paras,'SEM',1);  % SEM (for conventional mean) or STD (for bootstrap)
% addOptional(paras,'Border',[1600, -350]);
% 
% % addOptional(paras,'Markers',{'o','o','s','^','v','<','>'}); addOptional(paras,'MarkerSize',15);
% % addOptional(paras,'FaceColors',{'k','none',[0.8 0.8 0.8]});
% addOptional(paras,'Xlabel',[]);
% addOptional(paras,'Ylabel',[]);


% ------ Parse input parameters -------
paras = inputParser;

addOptional(paras,'OverrideError',[]); % Error overrided by inputs. HH20160908
addOptional(paras,'OverridePs',[]); % P values overrided by inputs. HH20160908

addOptional(paras,'figN',999);
addOptional(paras,'axes',[]);

addOptional(paras,'ErrorBar',2);  % Sum of (1: Normal Errorbar; 2: Shaded Errorbar; 4: Significance marker)
addOptional(paras,'CompareIndex',[]);  % Pair of statistic test for drawing significance marker: {[1 3 5; 2 4 6]}. HH20160427
addOptional(paras,'CompareColor',{});  % Colors for significance marker.
addOptional(paras,'PCritical',0.05);  % Colors for significance marker.

addOptional(paras,'Colors',{'b','r','g'}); % (originally for comparing PSTH of LIP neurons)
addOptional(paras,'LineStyles',{'-'});
addOptional(paras,'Transparent',1);
addOptional(paras,'SEM',1);  % SEM (for conventional mean) or STD (for bootstrap)
addOptional(paras,'Border',[1600, -350]);

% addOptional(paras,'Markers',{'o','o','s','^','v','<','>'}); addOptional(paras,'MarkerSize',15);
% addOptional(paras,'FaceColors',{'k','none',[0.8 0.8 0.8]});
addOptional(paras,'Xlabel',[]);
addOptional(paras,'Ylabel',[]);

parse(paras,varargin{:});

colors = paras.Results.Colors;
if ~iscell(colors)
    colors = mat2cell(colors,ones(size(colors,1),1));
end

compare_index = paras.Results.CompareIndex;
compare_color = paras.Results.CompareColor;
p_critical = paras.Results.PCritical;
border = paras.Results.Border;
override_error = paras.Results.OverrideError;
override_ps = paras.Results.OverridePs;

% --------- End input parser ----------

transparent = paras.Results.Transparent;

if isempty(paras.Results.axes)
    figure(paras.Results.figN); clf;
else
    axes(paras.Results.axes);
end
hold on;

if ~iscell(ys) % Only one j (old version)
    
    for cat = 1:size(ys,3)
        % Remove NaNs
        y = ys(:,:,cat);
        y(any(isnan(y),2),:) = [];
        
        % Plotting
        n = size(y,1);
        means = mean(y,1);
        
        % If overrided error
        if ~isempty(override_error)  % Overrided error. HH20160908
            errors = override_error(:,cat);
        else
            if paras.Results.SEM
                errors = std(y,[],1) / sqrt(n);
            else
                errors = std(y,[],1); % For bootstrapping
            end
        end
        
        error_type = find(fliplr(dec2bin(paras.Results.ErrorBar))=='1');
        
        if sum(error_type == 2)>0 && ~all(errors == 0)
            if ~all(isnan(errors))
                result.h(cat) = shadedErrorBar(ts,means,errors,...
                    {'Color',colors{1+mod(cat-1,length(colors))},...
                    'LineStyle',paras.Results.LineStyles{1+mod(cat-1,length(paras.Results.LineStyles))}},...
                    transparent);
                set(result.h(cat).mainLine,'LineWidth',2);
            end
        elseif sum(error_type == 1)>0
            result.h(cat) = errorbar(ts,means,errors,...
                'Color',colors{1+mod(cat-1,length(colors))},...
                'LineStyle',paras.Results.LineStyles{1+mod(cat-1,length(paras.Results.LineStyles))},'LineWidth',2.5);
        else 
            result.h(cat) = plot(ts,means,...
                'Color',colors{1+mod(cat-1,length(colors))},...
                'LineStyle',paras.Results.LineStyles{1+mod(cat-1,length(paras.Results.LineStyles))},'LineWidth',2.5);
        end
        
        % Saving data
        result.means(cat,:) = means;
        result.errors(cat,:) = errors;
        result.ns(cat) = n;
    end
    
    legend(num2str(result.ns'),'Location','best');
    xlabel(paras.Results.Xlabel);
    ylabel(paras.Results.Ylabel);
    if isempty(paras.Results.axes) ;SetFigure(); end
    
else % Deal with combined plot for two temporal alignments.  @HH20150523
    
    % Config
    %     border = [2500, -2500];
    gap = 50; % in ms
    time_markers = ts{3};
    
    % Plot range
    plot_range{1} = -inf <= ts{1} & ts{1} <= border(1);
    offset{1} = 0;
    plot_range{2} = border(2) <= ts{2} & ts{2} <= inf;
    offset{2} = border(1) + gap + -border(2);

    for j = 1:2
                
        for cat = 1:size(ys{j},3)
            % Remove NaNs
            y = ys{j}(:,:,cat);
            y(any(isnan(y),2),:) = [];
            
            % Plotting
            n = size(y,1);
            means = mean(y,1);
            % errors = std(y,[],1)/sqrt(n);
            
            % If overrided error
            if ~isempty(override_error)  % Overrided error
                errors = override_error{j}(:,cat);
            else
                
                if paras.Results.SEM
                    errors = std(y,[],1) / sqrt(n);
                else
                    errors = std(y,[],1); % For bootstrapping
                end
            end
            
            error_type = find(fliplr(dec2bin(paras.Results.ErrorBar))=='1');
            
            if sum(error_type == 2)>0 && (~all(errors == 0) && ~all(isnan(errors)))
                result.h{j}(cat) = shadedErrorBar(ts{j}(plot_range{j})+offset{j},means(plot_range{j}),errors(plot_range{j}),...
                    {'Color',colors{1+mod(cat-1,length(colors))},...
                    'LineStyle',paras.Results.LineStyles{1+mod(cat-1,length(paras.Results.LineStyles))}},...
                    transparent);
                set(result.h{j}(cat).mainLine,'LineWidth',2);
            else
                 result.h{j}(cat) = plot(ts{j}(plot_range{j})+offset{j},means(plot_range{j}),...
                    'Color',colors{1+mod(cat-1,length(colors))},...
                    'LineStyle',paras.Results.LineStyles{1+mod(cat-1,length(paras.Results.LineStyles))},'LineWidth',2.5);
            end
            
            % Saving data
            result.means{j}(cat,:) = means;
            result.errors{j}(cat,:) = errors;
            result.ns{j}(cat) = n;
            
        end
    end
      
    % Significance indicators. HH20160427
    axis tight; ylims = ylim;
    
    for j = 1:2
        if sum(error_type == 3)>0 && ~all(errors == 0)
            for cc = 1:size(compare_index,2)
                
                % Compute p values
                y1 = ys{j}(:,:,compare_index(1,cc));
                
                if compare_index(1,cc) ~= compare_index(2,cc)
                    y2 = ys{j}(:,:,compare_index(2,cc)); % Compare two series
                else
                    y2 = zeros(size(y1)); % Compare with zero
                end
                
                if ~isempty(override_ps) % HH20160908
                    ps = override_ps{j}(:,cc)';
                else
                    [~,ps] = ttest(y1,y2); % Should be paired t-test because in this version of SeriesComparison,
                    % we assume that each row (first dimension) of ys belongs to the same cell. HH20160427
                end
                dy = mean(y1,1)-mean(y2,1);
                
                % Plotting
                indicator_pos = zeros(size(dy)) + (cc-1)*ylims(2)/50;
                indicator_pos (dy > 0) = indicator_pos (dy > 0) + ylims(2)*1.05;
                
                lineColor = compare_color{cc};
                plot(ts{j}(plot_range{j} & ps<p_critical)+offset{j},indicator_pos(plot_range{j} & ps<p_critical),...
                    's','Color',lineColor,'MarkerFaceColor',lineColor,'MarkerSize',5);
                
            end
        end
    end
    
    % Time markers
    axis tight;
    marker_for_time_markers{1} = {'-','-','--'};
    marker_for_time_markers{2} = {'--','--','-'};
    
    for j = 1:2
        for tt = 1:3
            if  ts{j}(find(plot_range{j}==1,1)) <= time_markers{j}(1,tt) && time_markers{j}(1,tt) <= ts{j}(find(plot_range{j}==1,1,'last'))
                plot([1 1] * time_markers{j}(1,tt) + offset{j},ylim,'k','linestyle',marker_for_time_markers{j}{tt},'linew',2);
            end
        end
    end
    
    legend(num2str(result.ns{j}'),'Location','best');
    xlabel(paras.Results.Xlabel);
    ylabel(paras.Results.Ylabel);
    if isempty(paras.Results.axes) ;SetFigure(); end
    
    
end


