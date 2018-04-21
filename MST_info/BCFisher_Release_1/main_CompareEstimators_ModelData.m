%% This code loads synthetic data, estimates Fisher information 
%% with the direct bias corrected estimator, as well as 3 decoding-based 
%% methods, and then makes the main figures of the paper
%%
% Kanitscheider*, Coen-Cagli*, Kohn, Pouget. "Measuring Fisher information 
% accurately in correlataed neural populations". PLoS Comp Biol
% 2015

%% Load data
load ModelData

% resp: the spike count data generated from the model
%** size(resp) = [50 neurons  x  2 orientations  x  100,000 trials]
% FI_TRUE: the ground truth Fisher information (based on the analytical model tuning curves and covariances) 
% FI_TRUE_SHUF: the ground truth I_shuf
% FI_TRUE_DIAG: the ground truth I_diag 


%% parameters for information estimation
ORI = [-7 0]*pi/180; %rad
DORI = [1 2];
or_corr=[1 2];
ds = diff(ORI(or_corr));
NR=20; % how many cross-validation splits for early stopping
NP=200; % how many experiments
NT = [50 100 250 500 1000 2500 4000]; % how many trials per experiment
NTT=numel(NT);
N = size(resp,3);


%% ESTIMATE INFORMATION WITH DIFFERENT METHODS

%%% Analytical bias correction
FIBC = NaN(NTT,NP);
varFIBC = NaN(NTT,NP);
FINBC = NaN(NTT,NP); % naive, i.e. not bias corrected
FIBCSHUF = NaN(NTT,NP);
FINBCSHUF = NaN(NTT,NP);
FIBCDIAG = NaN(NTT,NP);
for t = 1:NTT
    NVAL = NT(t);
    for p=1:NP
        rng(p)
        indtmp = randperm(N); %subsample data to desired size
        D1 = squeeze(resp(:,or_corr(1),indtmp(1:NVAL)))';
        D2 = squeeze(resp(:,or_corr(2),indtmp(1:NVAL)))';
        [FIBC(t,p) varFIBC(t,p) FINBC(t,p)] = BCFisher(D1,D2,ds);
        [FIBCSHUF(t,p) FINBCSHUF(t,p)] = BCFisherShuf(D1,D2,ds);
        FIBCDIAG(t,p) = BCFisherDiag(D1,D2,ds);
    end
end

%%% DECODER WITH EARLY STOPPING
FITR = NaN(NTT,NP);
FIVAL = NaN(NTT,NP);
FISHUF = NaN(NTT,NP);
FISHUFTR = NaN(NTT,NP);
FIDIAG = NaN(NTT,NP);
fracTR=1/3;
fracTE=1/3;
for t = 1:NTT
    NVAL = NT(t);
    for p=1:NP
        rng(p)
        indtmp = randperm(N); %subsample data to desired size
        D1 = squeeze(resp(:,or_corr(1),indtmp(1:NVAL)))';
        D2 = squeeze(resp(:,or_corr(2),indtmp(1:NVAL)))';
        [FIVAL(t,p), FITR(t,p)] = EarlyStopping(D1,D2,ds,fracTR,fracTE,NR,1);
        [FISHUF(t,p), FISHUFTR(t,p),FIDIAG(t,p)] = EarlyStoppingShufDiag(D1,D2,ds,fracTR,fracTE,NR,1);
    end
end

%%% Variational Bayes Logistic Regression
FIVBTR = NaN(NTT,NP);
FIVBVAL = NaN(NTT,NP);
fracTRVB = 2/3;
for t = 1:NTT
    NVAL = NT(t);
    for p=1:NP
        rng(p)
        indtmp = randperm(N); %subsample data to desired size
        D1 = squeeze(resp(:,or_corr(1),indtmp(1:NVAL)))';
        D2 = squeeze(resp(:,or_corr(2),indtmp(1:NVAL)))';
        [FIVBVAL(t,p), FIVBTR(t,p)] = VBLogReg(D1,D2,ds,fracTRVB,NR,1);
    end
end

%%% LEAVE ONE OUT CV WITH RIDGE REGRESSION 
%***** Warning: this one takes very long time to run
FIRDGTR = NaN(NTT,NP);
FIRDGVAL = NaN(NTT,NP);
for t = 1:NTT
    NVAL = NT(t);
    for p=1:NP
        rng(p)
        indtmp = randperm(N); %subsample data to desired size
        D1 = squeeze(resp(:,or_corr(1),indtmp(1:NVAL)))';
        D2 = squeeze(resp(:,or_corr(2),indtmp(1:NVAL)))';       
        [FIRDGVAL(t,p), FIRDGTR(t,p)] = RidgeRegLOOCV(D1,D2,ds,0.1);
    end
end

%% Generate figures (Figure 1b-d; Figure 2a; Figure 4; Figure 5)

%*** I LOLE - BIAS CORRECTED vs DECODER WITH EARLY STOPPING
figure; 
subplot(1,2,1);hold on % INFORMATION ESTIMATE
myerrorbar(NT,nanmean(FINBC,2),nanstd(FINBC,[],2),[.85 .85 1],1)
myerrorbar(NT,nanmean(FIVAL,2),nanstd(FIVAL,[],2),[1 .85 .85],1)
myerrorbar(NT,nanmean(FITR,2),nanstd(FITR,[],2),[1 .85 .85],1)
myerrorbar(NT,nanmean(FIBC,2),nanstd(FIBC,[],2),[.65 .65 1],1)
plot(NT([1 end]),FI_TRUE([1 1]),'--k','LineWidth',2);
plot(NT,nanmean(FIVAL(:,1:NP),2),'-r')
plot(NT,nanmean(FITR(:,1:NP),2),'-r')
plot(NT,nanmean(FIBC(:,1:NP),2),'-b','LineWidth',2)
plot(NT,nanmean(FINBC(:,1:NP),2),'-b')
xlabel('Number of Trials')
ylabel('FI')
set(gca,'TickDir','out','XTick',NT,'xscale','log','XLim',[NT(1)*.75 NT(end)*1.5],'YLim',FI_TRUE*[0.25 2]); axis square
subplot(1,2,2); hold on % INFORMATION ESTIMATE ERROR
mnFIBC = nanmean(FIBC,2);
vrFIBC = nanvar(FIBC,[],2);
mnFI = nanmean(FIVAL,2);
vrFI = nanvar(FIVAL,[],2);
plot(NT,(mnFI-FI_TRUE).^2 + vrFI , '-r')
plot(NT,(mnFIBC-FI_TRUE).^2 + vrFIBC , '-b')
sqrt(10.^(1:5))/FI_TRUE %%TickLabels for percent error
xlabel('Number of Trials')
ylabel('Error(FI)')
set(gca,'TickDir','out','XTick',NT,'xscale','log','yscale','log','XLim',[NT(1)*.75 NT(end)*1.5]); axis square

%*** I SHUFFLE  - BIAS CORRECTED vs DECODER WITH EARLY STOPPING
figure; 
subplot(1,2,1);hold on % INFORMATION ESTIMATE
myerrorbar(NT,nanmean(FINBCSHUF,2),nanstd(FINBCSHUF,[],2),[.85 .85 1],1)
myerrorbar(NT,nanmean(FISHUF,2),nanstd(FISHUF,[],2),[1 .85 .85],1)
myerrorbar(NT,nanmean(FISHUFTR,2),nanstd(FISHUFTR,[],2),[1 .85 .85],1)
myerrorbar(NT,nanmean(FIBCSHUF,2),nanstd(FIBCSHUF,[],2),[.65 .65 1],1)
plot(NT([1 end]),FI_TRUE_SHUF*([1 1]),'--k','LineWidth',2);
plot(NT,nanmean(FISHUF,2),'-r')
plot(NT,nanmean(FISHUFTR,2),'-r')
plot(NT,nanmean(FIBCSHUF,2),'-b','LineWidth',2)
plot(NT,nanmean(FINBCSHUF,2),'-b')
xlabel('Number of Trials')
ylabel('FI')
set(gca,'TickDir','out','XTick',NT,'xscale','log','XLim',[NT(1)*.75 NT(end)*1.5],'YLim',FI_TRUE_SHUF*[0.25 2]); axis square
subplot(1,2,2); hold on % INFORMATION ESTIMATE ERROR
mnFIBC = nanmean(FIBCSHUF,2);
vrFIBC = nanvar(FIBCSHUF,[],2);
mnFI = nanmean(FISHUF,2);
vrFI = nanvar(FISHUF,[],2);
plot(NT,(mnFI-FI_TRUE_SHUF).^2 + vrFI , '-r')
plot(NT,(mnFIBC-FI_TRUE_SHUF).^2 + vrFIBC , '-b')
sqrt(10.^(1:5))/FI_TRUE_SHUF %%TickLabels for percent error
xlabel('Number of Trials')
ylabel('Error(FI)')
set(gca,'TickDir','out','XTick',NT,'xscale','log','yscale','log','XLim',[NT(1)*.75 NT(end)*1.5]); axis square


%*** I DIAG - BIAS CORRECTED vs DECODER WITH EARLY STOPPING
figure; 
subplot(1,2,1);hold on % INFORMATION ESTIMATE
myerrorbar(NT,nanmean(FIDIAG,2),nanstd(FIDIAG,[],2),[1 .85 .85],1)
myerrorbar(NT,nanmean(FIBCDIAG,2),nanstd(FIBCDIAG,[],2),[.65 .65 1],1)
plot(NT([1 end]),FI_TRUE_DIAG*([1 1]),'--k','LineWidth',2);
plot(NT,nanmean(FIDIAG,2),'-r')
plot(NT,nanmean(FIBCDIAG,2),'-b','LineWidth',2)
xlabel('Number of Trials')
ylabel('FI')
set(gca,'TickDir','out','XTick',NT,'xscale','log','XLim',[NT(1)*.75 NT(end)*1.5],'YLim',FI_TRUE_DIAG*[0.25 2]); axis square
subplot(1,2,2); hold on % INFORMATION ESTIMATE ERROR
mnFIBC = nanmean(FIBCDIAG,2);
vrFIBC = nanvar(FIBCDIAG,[],2);
mnFI = nanmean(FIDIAG,2);
vrFI = nanvar(FIDIAG,[],2);
plot(NT,(mnFI-FI_TRUE_DIAG).^2 + vrFI , '-r')
plot(NT,(mnFIBC-FI_TRUE_DIAG).^2 + vrFIBC , '-b')
sqrt(10.^(1:5))/FI_TRUE_DIAG %%TickLabels for percent error
xlabel('Number of Trials')
ylabel('Error(FI)')
set(gca,'TickDir','out','XTick',NT,'xscale','log','yscale','log','XLim',[NT(1)*.75 NT(end)*1.5]); axis square

%% Generate control Figure S7a,b

%*** I LOLE VAR BAYES
figure; 
subplot(1,2,1);hold on % INFORMATION ESTIMATE
myerrorbar(NT,nanmean(FIVBVAL,2),nanstd(FIVBVAL,[],2),[1 .85 .85],1)
myerrorbar(NT,nanmean(FIVBTR,2),nanstd(FIVBTR,[],2),[1 .85 .85],1)
myerrorbar(NT,nanmean(FIBC,2),nanstd(FIBC,[],2),[.65 .65 1],1)
plot(NT([1 end]),FI_TRUE*[1 1],'--k','LineWidth',2);
plot(NT,nanmean(FIVBVAL,2),'-r')
plot(NT,nanmean(FIVBTR,2),'-r')
plot(NT,nanmean(FIBC,2),'-b','LineWidth',2)
xlabel('Number of Trials')
ylabel('FI')
set(gca,'TickDir','out','XTick',NT,'xscale','log','XLim',[NT(1)*.75 NT(end)*1.5],'YLim',FI_TRUE*[0.25 2]); axis square
subplot(1,2,2); hold on % INFORMATION ESTIMATE ERROR
mnFIBC = nanmean(FIBC,2);
vrFIBC = nanvar(FIBC,[],2);
mnFI = nanmean(FIVBVAL,2);
vrFI = nanvar(FIVBVAL,[],2);
plot(NT,(mnFI-FI_TRUE).^2 + vrFI , '-r')
plot(NT,(mnFIBC-FI_TRUE).^2 + vrFIBC , '-b')
sqrt(10.^(1:5))/FI_TRUE %%TickLabels for percent error
xlabel('Number of Trials')
ylabel('Error(FI)')
set(gca,'TickDir','out','XTick',NT,'xscale','log','yscale','log','XLim',[NT(1)*.75 NT(end)*1.5]); axis square

%*** I LOLE RIDGE LOOCV
figure; 
subplot(1,2,1);hold on % INFORMATION ESTIMATE
myerrorbar(NT,nanmean(FIRDGVAL,2),nanstd(FIRDGVAL,[],2),[1 .85 .85],1)
myerrorbar(NT,nanmean(FIRDGTR,2),nanstd(FIRDGTR,[],2),[1 .85 .85],1)
myerrorbar(NT,nanmean(FIBC,2),nanstd(FIBC,[],2),[.65 .65 1],1)
plot(NT([1 end]),FI_TRUE*[1 1],'--k','LineWidth',2);
plot(NT,nanmean(FIRDGVAL,2),'-r')
plot(NT,nanmean(FIRDGTR,2),'-r')
plot(NT,nanmean(FIBC,2),'-b','LineWidth',2)
xlabel('Number of Trials')
ylabel('FI')
set(gca,'TickDir','out','XTick',NT,'xscale','log','XLim',[NT(1)*.75 NT(end)*1.5],'YLim',FI_TRUE*[0.25 2]); axis square
subplot(1,2,2); hold on % INFORMATION ESTIMATE ERROR
mnFIBC = nanmean(FIBC,2);
vrFIBC = nanvar(FIBC,[],2);
mnFI = nanmean(FIRDGVAL,2);
vrFI = nanvar(FIRDGVAL,[],2);
plot(NT,(mnFI-FI_TRUE).^2 + vrFI , '-r')
plot(NT,(mnFIBC-FI_TRUE).^2 + vrFIBC , '-b')
sqrt(10.^(1:5))/FI_TRUE %%TickLabels for percent error
xlabel('Number of Trials')
ylabel('Error(FI)')
set(gca,'TickDir','out','XTick',NT,'xscale','log','yscale','log','XLim',[NT(1)*.75 NT(end)*1.5]); axis square


%%


