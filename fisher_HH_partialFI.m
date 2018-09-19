function result = fisher_HH_partialFI(activities, headings, choices)

unique_heading = unique(headings);

infoPartialSensory = nan(size(activities,2),1);
infoPartialChoice = nan(size(activities,2),1);

% Fit cell by cell
for cc = 1:size(activities,2)
            
    % -  Multivariate linear regression of beta ---
    normalizedChoices = choices./rms(choices).*rms(headings); % Normalize choice to headings' RMS (to keep sensory heading comparable with the traditional FI)
    coeff = glmfit([headings normalizedChoices], activities(:,cc));
    
    conditional_sensory_slope = coeff(2) * (180/pi); % Turn to rad
    conditional_choice_slope = coeff(3) * (180/pi); % Turn to rad, although it's choice
    
    % --  Compute conditional variance
    conditional_variance_matrix = nan(2,length(unique_heading));
    unique_choice = unique(choices);
    
    for hh = 1:length(unique_heading)
        for ch = 1:length(unique_choice)
            this_rates = activities ( headings == unique_heading(hh) & choices == unique_choice(ch), cc);
            if length(this_rates) >= 5
                conditional_variance_matrix(ch,hh) = var(this_rates);
            end
        end
    end
    
    conditional_variance_mean = nanmean(conditional_variance_matrix(:));
    
    infoPartialSensory(cc) = conditional_sensory_slope^2/conditional_variance_mean;
    infoPartialChoice(cc) = conditional_choice_slope^2/conditional_variance_mean;
        
end

result.infoPartialSensory = nansum(infoPartialSensory);
result.bootSESensory = std(bootstrp(1000,@nansum,infoPartialSensory));
result.infoPartialChoice = nansum(infoPartialChoice);
result.bootSEChoice = std(bootstrp(1000,@nansum,infoPartialChoice));


