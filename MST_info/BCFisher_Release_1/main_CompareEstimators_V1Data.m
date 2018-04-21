%% This code loads the V1 data, estimates Fisher information 
%% with the direct bias corrected estimator, as well as with a decoder 
%% with early stopping, and then makes main Figure 6 of the paper
%%
% Kanitscheider*, Coen-Cagli*, Kohn, Pouget. "Measuring Fisher information 
% accurately in correlataed neural populations". PLoS Comp Biol
% 2015


%% load data and format
load V1data

[N K] = size(V1data_ORI0_NOISE0);
resp = NaN(N,3,2,K); % combine thespike counts from different conditions in one big matrix
resp(:,1,1,:) = V1data_ORI0_NOISE0;
resp(:,2,1,:) = V1data_ORI1_NOISE0;
resp(:,3,1,:) = V1data_ORI2_NOISE0;
resp(:,1,2,:) = V1data_ORI0_NOISE1;
resp(:,2,2,:) = V1data_ORI1_NOISE1;
resp(:,3,2,:) = V1data_ORI2_NOISE1;

ORI = [ORI0 ORI1 ORI2]; % stimulus orientations, in radians

%% parameters for information estimation
DORI = [[1 2]; [2 3]; [1 3]]; 
NDOR = size(DORI,1);
NR=20; % cross-validation splits
NP=100; % subsampling trials
NT = [50 75 100 150 300 600 900];
NTT=numel(NT);
NSIG=2; %levels of external noise
NPR = NR*NP;
subsample =1;

%% ESTIMATE INFORMATION WITH DIFFERENT METHODS

%%% Analytical bias correction
FIBC = NaN(NDOR,NSIG,NTT,NP);
varFIBC = NaN(NDOR,NSIG,NTT,NP);
FINBC = NaN(NDOR,NSIG,NTT,NP); % naive, i.e. not bias corrected
FIBCSHUF = NaN(NDOR,NSIG,NTT,NP);
FINBCSHUF = NaN(NDOR,NSIG,NTT,NP);
FIBCDIAG = NaN(NDOR,NSIG,NTT,NP);
for ndor = 1:NDOR
    or_corr=DORI(ndor,:);
    ds = diff(ORI(or_corr));
    for s = 1:NSIG
        for t = 1:NTT
            NVAL = NT(t);
            for p=1:NP
                rng(p)
                indtmp = randperm(N); %subsample data to desired size
                D1 = squeeze(resp(indtmp(1:NVAL),or_corr(1),s,:));
                D2 = squeeze(resp(indtmp(1:NVAL),or_corr(2),s,:));
                [FIBC(ndor,s,t,p) varFIBC(ndor,s,t,p) FINBC(ndor,s,t,p)] = BCFisher(D1,D2,ds);
                [FIBCSHUF(ndor,s,t,p) FINBCSHUF(ndor,s,t,p)] = BCFisherShuf(D1,D2,ds);
                FIBCDIAG(ndor,s,t,p) = BCFisherDiag(D1,D2,ds);
            end
        end
    end
end

%%% DECODER WITH EARLY STOPPING
FITR = NaN(NDOR,NSIG,NTT,NP);
FIVAL = NaN(NDOR,NSIG,NTT,NP);
FISHUF = NaN(NDOR,NSIG,NTT,NP);
FISHUFTR = NaN(NDOR,NSIG,NTT,NP);
FIDIAG = NaN(NDOR,NSIG,NTT,NP);
fracTR=1/3;
fracTE=1/3;
for ndor = 1:NDOR
    or_corr=DORI(ndor,:);
    ds = diff(ORI(or_corr));
    for s = 1:NSIG
        for t = 1:NTT
            NVAL = NT(t);
            for p=1:NP
                rng(p)
                indtmp = randperm(N); %subsample data to desired size
                D1 = squeeze(resp(indtmp(1:NVAL),or_corr(1),s,:));
                D2 = squeeze(resp(indtmp(1:NVAL),or_corr(2),s,:));
                [FIVAL(ndor,s,t,p), FITR(ndor,s,t,p)] = EarlyStopping(D1,D2,ds,fracTR,fracTE,NR,1);
                [FISHUF(ndor,s,t,p), FISHUFTR(ndor,s,t,p),FIDIAG(ndor,s,t,p)] = EarlyStoppingShufDiag(D1,D2,ds,fracTR,fracTE,NR,1);
            end
        end
    end
end

%% Generate figures (Figure 6, Figure S4)

%%% BC vs CV
nsamp=1000
figure; % information in full population, compare estimators
iplot=0;
for ndor=[1 NDOR]
    for s=[1:NSIG]
        iplot = iplot+1;
        subplot(2,2,iplot); hold on; axis square
        or_corr=DORI(ndor,:);
        myerrorbar(NT,squeeze(nanmean(FITR(ndor,s,:,:),4)),nanmean(bootstrp(nsamp,@nanstd,squeeze(FITR(ndor,s,:,:))')),[1 .85 .85],1);
        myerrorbar(NT,squeeze(nanmean(FIVAL(ndor,s,:,:),4)),nanmean(bootstrp(nsamp,@nanstd,squeeze(FIVAL(ndor,s,:,:))')),[1 .85 .85],1);
        myerrorbar(NT,squeeze(nanmean(FIBC(ndor,s,:,:),4)), nanmean(bootstrp(nsamp,@nanstd,squeeze(FIBC(ndor,s,:,:))')),[.85 .85 1],1);
        plot(NT,squeeze(nanmean(FIBC(ndor,s,:,:),4)),'-b','LineWidth',2);
        plot(NT,squeeze(nanmean(FITR(ndor,s,:,:),4)),'-r');
        plot(NT,squeeze(nanmean(FIVAL(ndor,s,:,:),4)),'-r');
        set(gca,'TickDir','out','xscale','log','XLim',[40 1000],'XTick',NT)
    end
end

%%% BC vs CV, mean squared error
nsamp=1000
figure; % squared error, compare estimators
iplot=0;
for ndor=[1 NDOR]
    for s=[1:NSIG]
        iplot = iplot+1;
        subplot(2,2,iplot); hold on; axis square
        or_corr=DORI(ndor,:);
        FIASYMPT = 0.5*(nanmean(FITR(ndor,s,end,:),4)+nanmean(FIVAL(ndor,s,end,:),4));
        
        mnFIBC = squeeze(nanmean(FIBC(ndor,s,:,:),4));
        vrFIBC = squeeze(nanvar(FIBC(ndor,s,:,:),[],4));
        mnFI = squeeze(nanmean(FIVAL(ndor,s,:,:),4));
        vrFI = squeeze(nanvar(FIVAL(ndor,s,:,:),[],4));
        plot(NT,(mnFI-FIASYMPT).^2 + vrFI , '-r')
        plot(NT,(mnFIBC-FIASYMPT).^2 + vrFIBC , '-b')
        sqrt(10.^(1:5))/FIASYMPT %%TickLabels for percent error
        xlabel('Number of Trials')
        ylabel('Error(FI)')
        set(gca,'TickDir','out','XTick',NT,'xscale','log','yscale','log','XLim',[NT(1)*.75 NT(end)*1.5]); axis square
    end
end


%% I shuf and I diag, Figure S5

% I shuf
figure; % information in shuffled population, compare estimators
iplot=0;
for ndor=[1 NDOR]
    for s=[1:NSIG]
        iplot = iplot+1;
        subplot(2,2,iplot); hold on; axis square
        or_corr=DORI(ndor,:);
        myerrorbar(NT,squeeze(nanmean(FIBCSHUF(ndor,s,:,:),4)), nanmean(bootstrp(nsamp,@nanstd,squeeze(FIBCSHUF(ndor,s,:,:))')),[.85 .85 1],1);
        myerrorbar(NT,squeeze(nanmean(FISHUFTR(ndor,s,:,:),4)),nanmean(bootstrp(nsamp,@nanstd,squeeze(FISHUFTR(ndor,s,:,:))')),[1 .85 .85],1);
        myerrorbar(NT,squeeze(nanmean(FISHUF(ndor,s,:,:),4)),nanmean(bootstrp(nsamp,@nanstd,squeeze(FISHUF(ndor,s,:,:))')),[1 .85 .85],1);
        plot(NT,squeeze(nanmean(FIBCSHUF(ndor,s,:,:),4)),'-b','LineWidth',2);
        plot(NT,squeeze(nanmean(FISHUFTR(ndor,s,:,:),4)),'-r');
        plot(NT,squeeze(nanmean(FISHUF(ndor,s,:,:),4)),'-r');
        set(gca,'TickDir','out','xscale','log','XLim',[40 1000],'XTick',NT)
    end
end

% I diag
figure; % information in shuffled population, compare estimators
iplot=0;
for ndor=[1 NDOR]
    for s=[1:NSIG]
        iplot = iplot+1;
        subplot(2,2,iplot); hold on; axis square
        or_corr=DORI(ndor,:);
        myerrorbar(NT,squeeze(nanmean(FIBCSHUF(ndor,s,:,:),4)), nanmean(bootstrp(nsamp,@nanstd,squeeze(FIBCSHUF(ndor,s,:,:))')),[.85 .85 1],1);
        myerrorbar(NT,squeeze(nanmean(FISHUFTR(ndor,s,:,:),4)),nanmean(bootstrp(nsamp,@nanstd,squeeze(FISHUFTR(ndor,s,:,:))')),[1 .85 .85],1);
        myerrorbar(NT,squeeze(nanmean(FISHUF(ndor,s,:,:),4)),nanmean(bootstrp(nsamp,@nanstd,squeeze(FISHUF(ndor,s,:,:))')),[1 .85 .85],1);
        plot(NT,squeeze(nanmean(FIBCSHUF(ndor,s,:,:),4)),'-b','LineWidth',2);
        plot(NT,squeeze(nanmean(FISHUFTR(ndor,s,:,:),4)),'-r');
        plot(NT,squeeze(nanmean(FISHUF(ndor,s,:,:),4)),'-r');
        set(gca,'TickDir','out','xscale','log','XLim',[40 1000],'YLim',[.25*min(squeeze(nanmean(FISHUF(ndor,s,:,:),4))) 1.25*max(squeeze(nanmean(FISHUFTR(ndor,s,:,:),4)))],'XTick',NT)
    end
end


%%