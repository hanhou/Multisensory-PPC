function [info, svm_model_final] = fisher_HH(activities, headings, train_ratio, svm_model_override)
%%  Use SVM decoder to calculate information in the population activity
%   [w,info] = fisher_HH(activities, headings, train_ratio)
%   activities: [trials, neurons]
%   headings: [trials, 1]
info_cutoff = 100;

if_debug = 0;
figN = round(rand()*100);

% Ignore zero headings (or randomly assign directions)
activities(headings==0,:) = [];
headings(headings==0) = [];

unique_heading = unique(headings);
n_rep = sum(headings==unique_heading(1)); % Assume we have balanced heading conditions

if nargin == 3 % No svm_model_override
    
    % To make sure each heading always has the same number of trials.
    n_train_each_condition = round(n_rep * train_ratio);
    n_test_each_condition = n_rep - n_train_each_condition;
    
    % Generate train and test sets
    train_set = false(size(headings));
    for hh = 1:length(unique_heading)
        train_this = randperm(n_rep,n_train_each_condition);
        train_set(train_this + (hh-1)*n_rep) = true;
    end
    
    X_train = activities(train_set,:);
    targ_train = sign(headings(train_set));
    
    X_test = activities(~train_set,:);
    targ_test = sign(headings(~train_set));
    
    % Train linear SVM, try different Cs, alike the "early stopping" method
    Cs = 10.^(-4:0.2:-1);
    % Cs = 1e-5;
    correct_rate_svm_test_Cs = nan(1,length(Cs));
    correct_rate_svm_train_Cs = correct_rate_svm_test_Cs;
    
    for ccs = 1:length(Cs)
        try
            svm_model_Cs{ccs} = svmtrain(X_train,targ_train,'boxconstraint',Cs(ccs),'tol',1e-7);
            
            choices_svm_train_Cs{ccs} = svmclassify(svm_model_Cs{ccs},X_train);
            correct_rate_svm_train_Cs(ccs) = sum(choices_svm_train_Cs{ccs} == targ_train) / length(targ_train);
            
            choices_svm_test_Cs{ccs} = svmclassify(svm_model_Cs{ccs},X_test);
            correct_rate_svm_test_Cs(ccs) = sum(choices_svm_test_Cs{ccs} == targ_test) / length(targ_test);
        catch
            % keyboard
        end
    end
    
    % Select the best C
    if if_debug
        figure(figN); clf; subplot(1,2,1);
        plot(correct_rate_svm_train_Cs,'o-k'); hold on; plot(correct_rate_svm_test_Cs,'o-r')
        legend('train','test');
    end
    
    bestC_ind = find(correct_rate_svm_test_Cs == max(correct_rate_svm_test_Cs));
    
    if length(bestC_ind)==1
        bestC = Cs(bestC_ind);
        svm_model_final = svm_model_Cs{bestC_ind};
        choices_svm_test_final = choices_svm_train_Cs{bestC_ind};
    else
        bestC = nanmean(Cs(bestC_ind));
        svm_model_final = svmtrain(X_train,targ_train,'boxconstraint',bestC,'tol',1e-7);
        choices_svm_test_final = svmclassify(svm_model_final,X_test);
    end
    
    if if_debug
        title(['C = ' num2str(bestC)]);
    end
    
    % Compute psychometric curves and Fisher information
    psychometric_test = [unique_heading sum(reshape(choices_svm_test_final == 1,[],length(unique_heading)),1)'/n_test_each_condition];
    [bias_test, threshold_test] = cum_gaussfit_max1(psychometric_test);
    
    % Use all data by simply combining train and test sets (should do bootstrapping here?) to increase N_rep
    choices_svm_all_final = svmclassify(svm_model_final,activities);
    psychometric_all = [unique_heading sum(reshape(choices_svm_all_final == 1,[],length(unique_heading)),1)'/n_rep];
    [bias_all, threshold_all] = cum_gaussfit_max1(psychometric_all);
    
    info = 2 * 1/threshold_test^2;  % Fisher info = (d'/ds)^2  %% <Theoretical Neuroscience Equ 3.49> 
                                    %             = (ds/sigma_r/ds)^2 = (1/sigma_r)^2 = (sqrt(2)/sigma_psycho)^2 = 2/threshold^2.
    if info > info_cutoff % || abs(bias_test) >= max(abs(headings))
        info = 0;
    end
    
    % info = 2 * 1/threshold_all^2;  % Fisher info = (d'/ds)^2 = (ds/sigma_r/ds)^2 = (1/sigma_r)^2 = (sqrt(2)/sigma_psycho)^2 = 2/threshold^2.
    
    if if_debug
        figure(figN); 
        subplot(1,2,2); hold on
        xxx = -max(unique_heading):0.1:max(unique_heading);
        
        plot(psychometric_test(:,1),psychometric_test(:,2),'or'); % Psychometric
        plot(psychometric_all(:,1),psychometric_all(:,2),'ob'); % Psychometric
        legend('test','all');
        
        plot(xxx,normcdf(xxx,bias_test,threshold_test),'-','color','r','linewid',4);
        plot(xxx,normcdf(xxx,bias_all,threshold_all),'-','color','b','linewid',4);
        
        axis([-8.5 8.5 0 1]);
        
        file_name = sprintf('./result/6p5_Info_debug_%g',figN);
        try
            export_fig('-painters','-nocrop','-m1.5',[file_name '.png']);
            saveas(gcf,[file_name '.fig'],'fig');
        end

    end
    
else % Use svm_model_override
    svm_model_final = svm_model_override;
    
    % Use all data
    choices_svm_all_final = svmclassify(svm_model_final,activities);
    psychometric_all = [unique_heading sum(reshape(choices_svm_all_final == 1,[],length(unique_heading)),1)'/n_rep];
    [bias_all, threshold_all] = cum_gaussfit_max1(psychometric_all);
        
    info = 2 * 1/threshold_all^2;  % Fisher info = (d'/ds)^2 = (ds/sigma_r/ds)^2 = (1/sigma_r)^2 = (sqrt(2)/sigma_psycho)^2 = 2/threshold^2.

    if info > info_cutoff % || abs(bias_all) >= max(abs(headings))
        info = 0;
        %{
        figure(2008); 
        hold on
        xxx = -max(unique_heading):0.1:max(unique_heading);
        
        plot(psychometric_all(:,1),psychometric_all(:,2),'ob'); % Psychometric
        plot(xxx,normcdf(xxx,bias_all,threshold_all),'-','color','b','linewid',4);
        
        axis([-8.5 8.5 0 1]);
        
        file_name = sprintf('./result/6p5_Info_debug_%g',2008);
        export_fig('-painters','-nocrop','-m1.5',[file_name '.png']);
        saveas(gcf,[file_name '.fig'],'fig');
        %}
    end

end