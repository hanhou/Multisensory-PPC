rep = 10;

times = nan(rep,24);

for del_cor = 0:23 % From 24 workers to 1 worker
	delete(gcp('nocreate'));
	parpool('node03',24-del_cor);
	
	for r = 1:rep
		tic
		lip_HH
		times(r,del_cor+1)=toc;
	end
	
	eval(sprintf('!stopworker -jobmanager node03 -name node03.local_worker%02g',24-del_cor)); 
	pause(10);
end
