% rep = 1;
% 
% for del_cor = 0:14 % From 29 workers to 1 worker
% 	delete(gcp('nocreate'));
% 	parpool('node03',29-del_cor*2);
% 	
% 	for r = 1:rep
% 		lip_HH;
% 		times(del_cor+1,1,r)=par_time;
% 		times(del_cor+1,2,r)=total_time
% 	end
% 	
% 	eval(sprintf('!stopworker -clean -name node03.local_worker%02g',29-del_cor*2)); 
%         eval(sprintf('!stopworker -clean -name node03.local_worker%02g',29-del_cor*2-1));
% 
% 	pause(10);
% end
rep = 10;
times = [];
n_cores = [55:-1:10 1];

for nnnn = 1:length(n_cores)
	delete(gcp('nocreate'));
	parpool('node03',n_cores(nnnn));
    
	for r = 1:rep
		lip_HH;
		times(nnnn,r)=par_time
		% times(nnnn,2,r)=total_time;
	end
	nanmean(times,2)
end

save times times n_cores