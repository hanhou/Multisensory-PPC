clear

global theta_pref
global target;
nu_unit=1008;
step=180/(nu_unit/4);
for k=1:nu_unit/4
    for l=1:4
       theta_pref((k-1)*4+l)= (-90+step*(k-1))/180*pi;
    end
end

load sharpen1\dataout895905;
data1 = data895;
data2 = data905;
cov1=cov(data895);


mean1 = mean(data1);
var1 = var(data1);
mean2 = mean(data2);
var2 = var(data2);

%%% fits circ norm through the mean
target = mean1';
param=[max(target),1.5,theta_pref*target/sum(target),min(target)];
a = fminsearch(@dist_square,param);
mean_act_fit1 =  circ_norm(a(1),a(2),theta_pref,a(3),a(4));



var0 = (var1+var2)/2;
der = mean2-mean1;

info = sum(der.^2./(var0+1e-30))
fprintf('Info :  %f\n',info);


clear data1 data2

load nosharpen1\dataout895905;
data1 = data895;
data2 = data905;
cov2=cov(data895);

mean1_2 = mean(data1);
var1_2 = var(data1);
mean2_2 = mean(data2);
var2_2 = var(data2);

%%% fits circ norm through the mean
target = mean1_2';
param=[max(target),1.5,theta_pref*target/sum(target),min(target)];
a = fminsearch(@dist_square,param);
mean_act_fit2 =  circ_norm(a(1),a(2),theta_pref,a(3),a(4));


var0_2 = (var1_2+var2_2)/2;
der_2 = mean2_2-mean1_2;

info_2 = sum(der_2.^2./(var0_2+1e-30))
fprintf('Info :  %f\n',info_2);


theta_pref=theta_pref/(pi)*180;

subplot(221)
plot(theta_pref,mean1_2,'b*')
hold on
plot(theta_pref,mean1,'r*')
hold off
axis tight;

subplot(222)
[Y1,I1]=max(mean_act_fit1);
[Y2,I2]=max(mean_act_fit2);
offset=I2-I1;
%mean_act_fit1=(mean_act_fit1-min(mean_act_fit1))/max(mean_act_fit1-min(mean_act_fit1));
%mean_act_fit2=(mean_act_fit2-min(mean_act_fit2))/max(mean_act_fit2-min(mean_act_fit2));
mean_act_fit2=mean_act_fit2(mod(offset:nu_unit-1+offset,1007)+1);
plot(theta_pref,mean_act_fit1,'r')
hold on
plot(theta_pref,mean_act_fit2,'b')
hold off
axis('tight');

%%% mean and variance are scaled to account for the fact that we had data over 500ms only
subplot(223)
plot(log10(mean1_2),log10(var1_2/2),'b.');
hold on
plot(log10(mean1),log10(var1/2),'r.');
hold off
axis([-1 2 -1 2]);

n=1;
for k=1:nu_unit
   if log10(mean1_2(k))>=0.1
      X(n)= log10(mean1_2(k));
      Y(n)= log10(var1_2(k));
      n=n+1;
   end
end
[P,S]=polyfit(X,Y,1);
fprintf('a=%4.2f   b=%4.2f\n\n',P(1),10^P(2));
hold on
%plot(log10(mean1_2),P(1)*log10(mean1_2)+P(2),'g')
hold off
n=1;
for k=1:nu_unit
   if log10(mean1_2(k))>=0.1
      X(n)= log10(mean1(k));
      Y(n)= log10(var1(k));
      n=n+1;
   end
end
[P,S]=polyfit(X,Y,1);
fprintf('a=%4.2f   b=%4.2f\n\n',P(1),10^P(2));
hold on
%plot(log10(mean1_2),P(1)*log10(mean1_2)+P(2),'c')
hold off

subplot(224)
plot([8 16 32 63 126 252 504 1008],[0.0711 0.0948 0.2317 0.3307 0.6418 1.3068 2.7718 5.24],'r.-');
hold on 
plot([8 16 32 63 126 252 504 1008],[0.0832 0.1030 0.2500 0.3193 0.5425 0.8370 1.1169 1.4704],'r.--');
plot([8 16 32 63 126 252 504 1008],[0.0115 0.0298 0.0485 0.0498 0.0938 0.1238 0.1622 0.2216],'g.--');
plot([8 16 32 63 126 252 504 1008],[0.0131 0.0241 0.0527 0.0939 0.2053 0.3439 0.6924 1.879],'g.-');
hold off
axis tight;

load nosharpen2/error
[aux1 aux2]=size(sampled);
sampled=sampled(2:aux2);
[aux1 aux2] = min(sampled);
sampled_shap=sampled(1:aux2-1);

clear sampled
load sharpen2/error
[aux1 aux2]=size(sampled);
sampled=sampled(2:aux2);
[aux1 aux2] = min(sampled);
sampled_shar=sampled(1:aux2-1);

plot(sampled_shap);
hold on
plot(sampled_shar,'r');
axis([0 2100 0 max(sampled_shar)]);
hold off

for j=1:1008
   aux(j)=cov1(1009-j,j);
end
%plot(aux);