%% 
addpath(genpath('/home/hh/Multisensory-PPC/util'));

rhos = linspace(0.001,10,20);
f0 = 50; % Normal range of MST firing. Hz
dts  = 10.^(linspace(-2,0,20))'/f0; % f*dt from 0:1

err_firing = nan(length(rhos),length(dts));
err_fano = err_firing;
err_max_corr = err_firing;
mean_lessThan90 = err_firing; 

for rr = 1:length(rhos)
    disp(rr);
    
    parfor tt = 1:length(dts)
        result = generateCorrelatedPoisson(rhos(rr),f0,dts(tt));
        err_firing(rr,tt) = result.err_firing_rate;
        err_fano(rr,tt) = result.err_fano;
        err_max_corr(rr,tt) = result.err_max_corr;
        mean_lessThan90(rr,tt) = result.mean_corr_less_than_90;
    end
end

%% Plotting

fdt_Beck2008 = 20 * 1e-3;

figure(1402); clf; set(gcf,'uni','norm','pos',[0.065        0.12       0.659       0.739]);
subplot(2,2,1)
xaxis = log10(dts * f0);
imagesc(xaxis, rhos, err_firing/f0*100); axis xy;
hold on; plot(log10(fdt_Beck2008)*ones(1,2),ylim,'w--','linew',2);
contour(xaxis,rhos, err_firing/f0*100,'color','w','linew',1,'ShowText','on');
xlabel('log_{10}(\deltatf_0)')
ylabel('\rho')
title('(Actual-Desired)/Desired firing %');
colorbar;

 
subplot(2,2,2)
xaxis = log10(dts * f0);
imagesc(xaxis, rhos, err_fano); axis xy;
hold on; plot(log10(fdt_Beck2008)*ones(1,2),ylim,'w--','linew',2);
contour(xaxis,rhos, err_fano,'color','w','linew',1,'ShowText','on');
xlabel('log_{10}(\deltatf_0)')
ylabel('\rho')
title('Actual fano - 1');
colorbar;

subplot(2,2,3)
xaxis = log10(dts * f0);
imagesc(xaxis, rhos, err_max_corr); axis xy;
hold on; plot(log10(fdt_Beck2008)*ones(1,2),ylim,'w--','linew',2);
contour(xaxis,rhos, err_max_corr,'color','w','linew',1,'ShowText','on');
xlabel('log_{10}(\deltatf_0)')
ylabel('\rho')
title('Actual - Desired \rho');
colorbar;
colormap(jet)

subplot(2,2,4)
xaxis = log10(dts * f0);
imagesc(xaxis, rhos, mean_lessThan90); axis xy;
hold on; plot(log10(fdt_Beck2008)*ones(1,2),ylim,'w--','linew',2);
contour(xaxis,rhos, mean_lessThan90,'color','w','linew',1,'ShowText','on');
xlabel('log_{10}(\deltatf_0)')
ylabel('\rho')
title('Actual cor90');
colorbar;
colormap(jet)


saveas(gcf,'ScanCorrPoisson.fig');
