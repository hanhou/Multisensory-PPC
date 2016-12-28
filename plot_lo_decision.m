
clear

load performance
coherence=[1 3.2 6.4 12.8 25.6 51.2];
    
load diff_ac_cond1
diff_ac_mat(1,:)=diff_ac_cond;

load diff_ac_cond2
diff_ac_mat(2,:)=diff_ac_cond;

load diff_ac_cond3
diff_ac_mat(3,:)=diff_ac_cond;

load diff_ac_cond4
diff_ac_mat(4,:)=diff_ac_cond;

load diff_ac_cond5
diff_ac_mat(5,:)=diff_ac_cond;

load diff_ac_cond6
diff_ac_mat(6,:)=diff_ac_cond;

  
%errorbar(coherence,mean(diff_ac_mat),-std(diff_ac_mat),std(diff_ac_mat) );
semilogx(coherence,mean(diff_ac_mat),'s-');

%plot(100*performance,mean(diff_ac_mat),'s-')
%axis([50 100 0.9*min(mean(diff_ac_mat)) 1.1*max(mean(diff_ac_mat))])