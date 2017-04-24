

subplot(221)
wind_times = wind_offset+[0:nu_window]*wind_size;

if nu_targets==2
   lik1 = [1 pdf_aver(1,50) pdf_aver(2,50) pdf_aver(3,50) pdf_aver(4,50)];
   lik2 = [1 pdf_aver(1,100) pdf_aver(2,100) pdf_aver(3,100) pdf_aver(4,100)];
   loglik= log10(lik1./lik2);
%    if shuffle_flag ==1
%       lik1_shuff = [1 pdf_shuff_aver{1}(50) pdf_shuff_aver{2}(50) pdf_shuff_aver{3}(50) pdf_shuff_aver{4}(50)];
%       lik2_shuff = [1 pdf11_shuff_aver{1}(100) pdf_shuff_aver{2}(100) pdf_shuff_aver{4}(100) pdf_shuff_aver{4}(100)];
%       loglik_shuff= log10(lik1_shuff./lik2_shuff);
%    end
   
% size   plot(wind_times, loglik,'*-');
    plot(wind_times, loglik,'*-');
%    if shuffle_flag==1
%        hold on
%        plot([1 100 200 300 trial_dur], loglik_shuff,'r*-');
%        hold off
%    end
   title('Log odds');
elseif info_out_flag==1 
   plot(wind_times,[0 info_outx],'r-o')
   hold on
   if shuffle_flag==1
      plot(wind_times, [0 info_shuff_out],'b-o')
   end
   plot(wind_times,[0 info_train_outx],'r:o');
   plot(wind_times,[0 (info_outx+info_train_outx)/2],'g:o');
   plot(wind_times,ML_info,'c:*');
   hold off
   title('Fisher information');
end


%%% plot information over time in input layer
%    plot(wind_times,[0 info_inx],'r-o')
%    hold on
%    plot(wind_times,[0 info_train_inx],'r:o');
%    plot(wind_times,[0 (info_inx+info_train_inx)/2],'g:o');
%    hold off
%    title('Fisher information in input layer');
% 





subplot(222)
plot(pdf{1}(:,1),'.r');
hold on 
plot(pdf{2}(:,1),'.b');
plot(pdf{3}(:,1),'.g');
plot(pdf{4}(:,1),'.c');
plot(pdf_aver(1,:),'r')
plot(pdf_aver(2,:),'b')
plot(pdf_aver(3,:),'g')
plot(pdf_aver(4,:),'c')
hold off

axis tight
ylabel('P(s|r)')
xlabel('Stimulus')





subplot(223)



for m=1:nu_window
plot(pos_in,mean_in,'o')
hold on
plot(pos_out,mean_out,'or')
plot(pos_in,fit_in,'-b')
plot(pos_in,fit_out,'-r')
hold off
max_plot=max(max(mean_out),max(mean_in));
axis([0 180 0 max_plot*1.2]);
    aux11{m} = mean(proba_count{1,m}');
    aux11{m} = aux11{m}.*(aux11{m}>0);
%     aux11{m} = aux11{m}-min(aux11{m});
%     aux11{m} = aux11{m}/max(aux11{m});
end

plot(aux11{1},'or')
hold on
plot(aux11{2},'ob')
plot(aux11{3},'og')
plot(aux11{4},'oc')
hold off




subplot(224)
%surf([0:180/(N_out-1):180], [0:180/(N_out-1):180],cov_out-diag(diag(cov_out)))
%shading interp
aux_proba= [proba_out_av0__pretrial(:,pre_trial_dur-49:10:pre_trial_dur) proba_out_av1(:,1:10:trial_dur)];
%surfl(aux_proba)

plot([-49:10:trial_dur],aux_proba(50,:))
hold on
plot([-49:10:trial_dur],aux_proba(100,:),'g')
plot([-49:10:trial_dur],aux_proba(25,:),'c')
plot([-49:10:trial_dur],aux_proba(50,:)-(aux_proba(100,:).*(aux_proba(100,:)>0)),'r')
hold off
axis([-49 trial_dur 0 60] )
title('Conditioned on neuron')
 

%%% compute and plot Tin and Tout trajectories conditioned on response
 proba_out50_cor = proba_out50.*(perf_trial'*ones(1,trial_dur));
 proba_out100_cor = proba_out100.*(perf_trial'*ones(1,trial_dur));
 
 neuron_traj = sum(proba_out50_cor)/(sum(perf_trial));
 antineuron_traj = sum(proba_out100_cor)/(sum(perf_trial));
%  
 plot(neuron_traj)
 hold on
 plot(mean(proba_out50),'b--')
 plot(antineuron_traj,'g');
 plot(mean(proba_out100),'g--')
 hold off
axis([0 200 20 60])
ylabel('Activity (Hz)');
xlabel('Time (ms)');
title('Conditioned on correct')


%save exp_data spike_time50 spike_time100 perf_trial coherence RT
save datalast RT conf_level perf_trial perf mean_RT mean_RT_incor coherence neuron_traj antineuron_traj

beep
toc
