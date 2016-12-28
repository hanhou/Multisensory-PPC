clear

load nosharpen1\dataout895905;

data1=data895;
data2=data905;
cov1=cov(data1);
cov2=cov(data2);
[nu_trial nu_unit] = size(data1);

nu_column = nu_unit*2 + (nu_unit*(nu_unit+1))/2;
qdata1 = zeros(nu_unit,nu_column);
qdata2 = zeros(nu_unit,nu_column);

fprintf('Done loading data');

for j=1:nu_unit
    qdata1(j,1:nu_unit)=data1(j,1:nu_unit);
    qdata1(j,nu_unit+1:2*nu_unit)=data1(j,1:nu_unit).^2;
    qdata2(j,1:nu_unit)=data2(j,1:nu_unit);
    qdata2(j,nu_unit+1:2*nu_unit)=data2(j,1:nu_unit).^2;
end

count=2*nu_unit+1
for j=2:nu_unit
    qdata1(j,count:count+j-2)=cov1(j,1:j-1);
    qdata2(j,count:count+j-2)=cov2(j,1:j-1);
    if mod(j,100)==0
         j
    end
    count=count+j-1
end


Qdata895=qdata1;
Qdata905=qdata2;
save nosharpen1\Qdataout895905 Qdata895 Qdata905;