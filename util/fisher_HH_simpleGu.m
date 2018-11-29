function [infoSimpleGu, bootSE] = fisher_HH_simpleGu(activities, headings)


unique_heading = unique(headings);

% Get tunings
tuningMean = nan(length(unique_heading),size(activities,2));
tuningVar = tuningMean;

for hh = 1:length(unique_heading)
    thisAct = activities(headings == unique_heading(hh),:);
    tuningMean(hh,:) = mean(thisAct,1);
    tuningVar(hh,:) = var(thisAct,1);
end

% Fit cell by cell
infoEachCell = nan(size(activities,2),1);

for cc = 1:size(activities,2)
    linearFit = polyfit(unique_heading,tuningMean(:,cc),1);
    slopeInRad = linearFit(1)*(180/pi);
    varReal = mean(tuningVar(:,cc)); % Using real variance (usually much large)
    
    % Compute simple Fisher (simple method like Gu)
    infoEachCell(cc) = slopeInRad^2/varReal;
end

infoSimpleGu = nansum(infoEachCell);
bootSE = std(bootstrp(1000,@nansum,infoEachCell));
