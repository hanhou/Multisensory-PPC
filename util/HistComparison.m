function h = HistComparison(raw,varargin) % (figN,logx,logy,seperate)
% Histogram and data comparison. (originally for plotting Microstimulation data)
% 
%  h = HistComparison(raw,varargin)
% 
%     addOptional(paras,'figN',2499);
%     addOptional(paras,'MeanType','Mean');
% 
%     addOptional(paras,'CombinedIndex',[]);
% 
%     addOptional(paras,'Style','grouped');
% 
%     addOptional(paras,'EdgeColors',{'k','none',[0.8 0.8 0.8]});
%     addOptional(paras,'FaceColors',{'k','none',[0.8 0.8 0.8]});
%     addOptional(paras,'Xlabel','X');
% 
% @HH20150113

% ------ Parse input parameters -------
paras = inputParser;
addOptional(paras,'XCenters',[]);
addOptional(paras,'XNumbers',15);
addOptional(paras,'FigN',2699);

addOptional(paras,'MeanCombinedIndex',[]);
addOptional(paras,'MeanType','Mean');

addOptional(paras,'TTest',[]);

addOptional(paras,'Style','stacked');

addOptional(paras,'EdgeColors',{'k','k','k'});
addOptional(paras,'FaceColors',{'k','none',[0.8 0.8 0.8]});
addOptional(paras,'Xlabel','X');

addOptional(paras,'Axes',[]);

parse(paras,varargin{:});
% --------- End input parser ----------

%%
% Reoranize
if ~iscell(raw)
    raw = mat2cell(raw,size(raw,1),ones(1,size(raw,2)));
end

% Histgram
h.x_centers = paras.Results.XCenters;
if isempty(h.x_centers)
    all_min = min(cellfun(@min,raw));
    all_max = max(cellfun(@max,raw));
    h.x_centers = linspace(all_min,all_max,paras.Results.XNumbers);
end

TTest = paras.Results.TTest;

% Minor adjustment of xcenters
if ~isempty(TTest)
    [~,ii] = min(abs(TTest - (h.x_centers(1:end-1)+h.x_centers(2:end))/2));
    h.x_centers = h.x_centers + TTest - (h.x_centers(ii)+h.x_centers(ii+1))/2;
end

hist_counts = cell2mat(cellfun(@(x)hist(x,h.x_centers),raw,'Uniform',0)')';

% Plots
if isempty(paras.Results.Axes)
    set(figure(paras.Results.FigN),'position',[704 88 626 425]); clf;   hold on;
    h.ax = gca;
else
    h.ax = paras.Results.Axes;
    axes(h.ax); hold on;
end
    
h.bar = bar(h.x_centers',hist_counts,1,paras.Results.Style);

for i = 1:length(raw) % Raw data
    set(h.bar(i),'EdgeColor',paras.Results.EdgeColors{1+mod(i-1,length(paras.Results.EdgeColors))},...
        'FaceColor',paras.Results.FaceColors{1+mod(i-1,length(paras.Results.FaceColors))});
end

% Mean annotation
ylims = ylim;

% h.means
if strcmpi(paras.Results.Style,'grouped')
    if strcmpi(paras.Results.MeanType,'Mean')
        h.means = cellfun(@mean,raw);
    else strcmpi(paras.Results.MeanType,'Median')
        h.means = cellfun(@median,raw);
    end
    
    for i = 1:length(raw) % Raw data
        h.mean(i) = plot(h.means(i),ylims(2)*1.1,[paras.Results.EdgeColors{1+mod(i-1,length(paras.Results.EdgeColors))} 'v'],...
        'MarkerFaceColor',paras.Results.FaceColors{1+mod(i-1,length(paras.Results.FaceColors))},'markersize',14,'linew',2);
    end
    
elseif strcmpi(paras.Results.Style,'stacked')
    if strcmpi(paras.Results.MeanType,'Mean')
%         try
%             h.means =  mean([raw{:}]);
%         catch
            h.means =  mean(cell2mat(raw(:)));
%         end
    else strcmpi(paras.Results.MeanType,'Median')
%         try
%             h.means =  median([raw{:}]);
%         catch
            h.means =  median(cell2mat(raw(:)));
%         end
    end
    
    h.mean = plot(h.means,ylims(2)*1.05,[paras.Results.EdgeColors{1+mod(i-1,length(paras.Results.EdgeColors))} 'v'],...
        'MarkerFaceColor',paras.Results.EdgeColors{1+mod(i-1,length(paras.Results.EdgeColors))},'markersize',14,'linew',2);
end

if ~isempty(TTest)
    plot([TTest TTest],ylim,'k--');
    
    if strcmpi(paras.Results.Style,'grouped')
        [~,h.ps] = cellfun(@(x)ttest(x,TTest),raw);
    elseif strcmpi(paras.Results.Style,'stacked')
        [~,h.ps] = ttest(cell2mat(raw(:)),TTest);
    end
end

% if length(h.mean)>1
%     if ~isempty(TTest)
%         legend(h.mean,num2str([h.means; h.ps]'));
%     else
%         legend(h.mean,num2str([h.means]'));
%     end
% else
for i = 1:length(h.mean)
    if ~isempty(TTest)
        text(h.means(i),ylims(2) * 1.05,sprintf('\\mu = %3.3g\np = %3.3g',h.means(i),h.ps(i)),'color',paras.Results.EdgeColors{1+mod(i-1,length(paras.Results.EdgeColors))});
    else
        text(h.means(i),ylims(2) * 1.05,sprintf('\\mu = %3.3g',h.means(i)),'color',paras.Results.EdgeColors{1+mod(i-1,length(paras.Results.EdgeColors))});
    end
end

% SetFigure(15); axis tight;

    


