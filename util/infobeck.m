function[Itrain,Itest,w,pctrain,pctest]=info(train,test,dtheta,dt)

[ns,n]=size(train);  % assumes that ns is number of samples and n is number of neurons and size(train) == size(test);
                    % On input dt should be less than one to ensure numerical stability but should 
                    % also be as small as possible (recommended dt=0.1).  (decreasing dt increases run time)
                    
trbar = mean(train)';
cxxtr = cov(train);
tr2=train(ns/2+1:ns,:);
tr1=train(1:ns/2,:);
clear train

tebar = mean(test)';
cxxte = cov(test) + tebar*tebar' + trbar*trbar' - tebar*trbar' - trbar*tebar';
te2=test(ns/2+1:ns,:);
te1=test(1:ns/2,:);
clear test
% this cov is different from the first cov since the mean which is removed from the 
% test set is the mean computed from the training set.

% compute  covarnance of training set
trbar2=mean(tr2)';
trbar1=mean(tr1)';
cxytr = (trbar2-trbar1);

% compute covariances for test set.
tebar2=mean(te2)';
tebar1=mean(te1)';
cxyte = (tebar2-tebar1);

dt=dt/n  % this should yeild numerical stability of euler time stepping.
%dt=dt/max(eig(cxxtr));  % wthis *will* yeild numerical stability of euler time stepping.
k=0;
w=zeros(n,1);
dwdt=zeros(n,1);
dEtedt=-1;  % the derivative of the error function on the test set with respect to "time"
while(dEtedt<0)
   w=w+dt*dwdt;  % iterates to solve w(t) which will minimize training set error.
   dwdt = (cxytr-cxxtr*w);
   dEtedt = -(cxyte-cxxte*w)'*dwdt;   % computes exactly the derivative of the test set error  
                                      % and stops which this derivative is positive
   k=k+1;
end
iterations = k

bptrain = w'*(trbar2-trbar1)/dtheta;
vartrain = w'*(cov(tr1)+cov(tr2))*w/2;

bptest = w'*(tebar2-tebar1)/dtheta;
vartest = w'*(cov(te1)+cov(te2))*w/2;

Itrain = (bptrain)^2/(vartrain)  % only place where units matter is on information
Itest  = (bptest)^2/(vartest)
Ibar = (Itrain+Itest)/2
w=(bptest+bptrain)/(vartest+vartrain)*w;

pctrain = (sum(tr1*w<w'*trbar)+sum(tr2*w>w'*trbar))/ns
pctest  = (sum(te1*w<w'*trbar)+sum(te2*w>w'*trbar))/ns

dprime = sqrt(Ibar)*dtheta  % unitless d'
pccheck = 1/2+1/2*erf(dprime/2/sqrt(2))  

% dp = sqrt(Itrain)*dtheta;  % unitless d'
% pcchecktrain = 1/2+1/2*erf(dp/2/sqrt(2))  
% dp = sqrt(Itest)*dtheta;  % unitless d'
% pcchecktest = 1/2+1/2*erf(dp/2/sqrt(2))  
