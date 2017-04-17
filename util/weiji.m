%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% reorganize data for wei ji
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear



load exp_data0-3
[aux nu_trial]=size(spike_time100)
%nu_trial = floor(nu_trial/2);
spike_time = cell(1,nu_trial*12);

for j=1:nu_trial
   spike_time{1,j}= spike_time50{j};
   spike_time{1,nu_trial+j}= spike_time100{j};
end

A(1:nu_trial,1) = ones(nu_trial,1);
A(nu_trial+1:2*nu_trial,1) = zeros(nu_trial,1);
A(1:nu_trial,2)=perf_trial;
A(nu_trial+1:2*nu_trial,2)=perf_trial;
A(1:nu_trial,3)=RT;
A(nu_trial+1:2*nu_trial,3)=RT;
A(1:2*nu_trial,4)=ones(nu_trial*2,1)*coherence;


clear spike_time50 spike_time100 perf_trial RT
load exp_data32-3
for j=1:nu_trial
   spike_time{1,nu_trial*2+j}= spike_time50{j};
   spike_time{1,nu_trial*3+j}= spike_time100{j};
end
A(nu_trial*2+1:nu_trial*3,1) = ones(nu_trial,1);
A(nu_trial*3+1:nu_trial*4,1) = zeros(nu_trial,1);
A(nu_trial*2+1:nu_trial*3,2) = perf_trial;
A(nu_trial*3+1:nu_trial*4,2) = perf_trial;
A(nu_trial*2+1:nu_trial*3,3) = RT;
A(nu_trial*3+1:nu_trial*4,3) = RT;
A(nu_trial*2+1:nu_trial*4,4) = ones(nu_trial*2,1)*coherence;



clear spike_time50 spike_time100 perf_trial RT
load exp_data64-3
for j=1:nu_trial
   spike_time{1,nu_trial*4+j}= spike_time50{j};
   spike_time{1,nu_trial*5+j}= spike_time100{j};
end
A(nu_trial*4+1:nu_trial*5,1) = ones(nu_trial,1);
A(nu_trial*5+1:nu_trial*6,1) = zeros(nu_trial,1);
A(nu_trial*4+1:nu_trial*5,2) = perf_trial;
A(nu_trial*5+1:nu_trial*6,2) = perf_trial;
A(nu_trial*4+1:nu_trial*5,3) = RT;
A(nu_trial*5+1:nu_trial*6,3) = RT;
A(nu_trial*4+1:nu_trial*6,4) = ones(nu_trial*2,1)*coherence;


clear spike_time50 spike_time100 perf_trial RT
load exp_data128-3
for j=1:nu_trial
   spike_time{1,nu_trial*6+j}= spike_time50{j};
   spike_time{1,nu_trial*7+j}= spike_time100{j};
end
A(nu_trial*6+1:nu_trial*7,1) = ones(nu_trial,1);
A(nu_trial*7+1:nu_trial*8,1) = zeros(nu_trial,1);
A(nu_trial*6+1:nu_trial*7,2) = perf_trial;
A(nu_trial*7+1:nu_trial*8,2) = perf_trial;
A(nu_trial*6+1:nu_trial*7,3) = RT;
A(nu_trial*7+1:nu_trial*8,3) = RT;
A(nu_trial*6+1:nu_trial*8,4) = ones(nu_trial*2,1)*coherence;
save all_exp_data spike_time A

clear spike_time50 spike_time100 perf_trial RT
load exp_data256-3
for j=1:nu_trial
   spike_time{1,nu_trial*8+j}= spike_time50{j};
   spike_time{1,nu_trial*9+j}= spike_time100{j};
end
A(nu_trial*8+1:nu_trial*9,1) = ones(nu_trial,1);
A(nu_trial*9+1:nu_trial*10,1) = zeros(nu_trial,1);
A(nu_trial*8+1:nu_trial*9,2) = perf_trial;
A(nu_trial*9+1:nu_trial*10,2) = perf_trial;
A(nu_trial*8+1:nu_trial*9,3) = RT;
A(nu_trial*9+1:nu_trial*10,3) = RT;
A(nu_trial*8+1:nu_trial*10,4) = ones(nu_trial*2,1)*coherence;

clear spike_time50 spike_time100 perf_trial RT
load exp_data512-3
for j=1:nu_trial
   spike_time{1,nu_trial*10+j}= spike_time50{j};
   spike_time{1,nu_trial*11+j}= spike_time100{j};
end
A(nu_trial*10+1:nu_trial*11,1) = ones(nu_trial,1);
A(nu_trial*11+1:nu_trial*12,1) = zeros(nu_trial,1);
A(nu_trial*10+1:nu_trial*11,2) = perf_trial;
A(nu_trial*11+1:nu_trial*12,2) = perf_trial;
A(nu_trial*10+1:nu_trial*11,3) = RT;
A(nu_trial*11+1:nu_trial*12,3) = RT;
A(nu_trial*10+1:nu_trial*12,4) = ones(nu_trial*2,1)*coherence;

save -v6 all_exp_data2 spike_time A
