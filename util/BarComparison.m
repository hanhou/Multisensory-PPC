function result = BarComparison(data,varargin)
% Draw bar comparison, recursive (paired) t-test, etc.
% Input: data   column for varibles, row for observations
%
% begin at 20160607 HH
% 
% addOptional(paras,'figN',555);
% addOptional(paras,'axes',[]);
% 
% addOptional(paras,'Colors',{'b','r','g'}); % (originally for comparing PSTH of LIP neurons)
% addOptional(paras,'LineStyles',{'ko-'});
% addOptional(paras,'SEM',1);
% addOptional(paras,'PairTTest',1);
% addOptional(paras,'EqualVar',1);

% ------ Parse input parameters -------
paras = inputParser;

addOptional(paras,'figN',555);
addOptional(paras,'axes',[]);

addOptional(paras,'Colors',{'b','r','g'}); % (originally for comparing PSTH of LIP neurons)
addOptional(paras,'LineStyles',{'ko-'});
addOptional(paras,'SEM',1);
addOptional(paras,'PairTTest',1);
addOptional(paras,'EqualVar',1);

parse(paras,varargin{:});

Colors = paras.Results.Colors;
LineStyle = paras.Results.LineStyles;
sem = paras.Results.SEM;
user_axes = paras.Results.axes;
pair_ttest = paras.Results.PairTTest;
equal_var = paras.Results.EqualVar;

% --------- End input parser ----------


% Basic statistics
result.means = nanmean(data,1);
result.stds = nanstd(data,[],1);
result.sems = result.stds./sqrt(sum(~isnan(data),1));
n_var = size(data,2);

% Draw bars
if isempty(user_axes)
    figure(paras.Results.figN); clf; 
else
    axes(user_axes);
end
hold on; xlim([0.5 n_var+0.5]);

for i = 1:n_var
    col_ind = 1+mod(i-1,length(Colors));
    
    bar(i,result.means(i),0.7,'facecol',Colors{col_ind},'edgecol','none'); 
    if sem
        herr = errorbar(i,result.means(i),result.sems(i),'linestyle','none','linewidth',3,'col',Colors{col_ind});
    else
        herr = errorbar(i,result.means(i),result.stds(i),'linestyle','none','linewidth',3,'col',Colors{col_ind});
    end
    errorbar_tick(herr,13);
end
        
% Overdraw raw data
if ~ strcmp(LineStyle,'none')
    plot(1:n_var,data, LineStyle);
end

% Recursive paired t-test
ps_ttest = nan(size(data,2));
ps_vartest = nan(size(data,2));
name = 1:size(data,2);

for g1 = 1:size(data,2)
    for g2 = g1 + 1:size(data,2)
        % Paired ttest
        if pair_ttest
            [~,p] = ttest(data(:,g1),data(:,g2));
        else
            if equal_var
                [~,p] = ttest2(data(:,g1),data(:,g2));
            else
                [~,p] = ttest2(data(:,g1),data(:,g2),'Vartype','unequal');
            end
        end
        
        ps_ttest(g1,g2) = p;
        ps_ttest(g2,g1) = p;
        
        % Var test
        [~,p_var] = vartest2(data(:,g1),data(:,g2));
        ps_vartest(g1,g2) = p_var;
        ps_vartest(g2,g1) = p_var;
    end
end
%  sum(sum(ps<0.05))

result.ps_ttest = [[nan;name'] [name;ps_ttest] ];
result.ps_vartest = [[nan;name'] [name;ps_vartest] ];



