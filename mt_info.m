clear

N_in = 100;
pos_in=[0:180/N_in:179.9999];

r_spont_MT = 20;
b_pref = 0.4;
b_null = -0.2;
K_in = 4;

coherence = [1:100];

    
    
for j=1:100
   max_rate_in = r_spont_MT + b_pref*coherence(j)*10;
   b_in = r_spont_MT + b_null*coherence(j);
   proba_in(:,j) = (max_rate_in-b_in)*exp(K_in*(cos(pos_in/180*2*pi)-1))+b_in; 
   proba_in_der(:,j) =(max_rate_in-b_in)*exp(K_in*(cos(pos_in/180*2*pi)-1))*...
       K_in.*-sin(pos_in/180*2*pi)/180*2*pi; 

end

info = sum(proba_in_der.^2./proba_in);

plot(coherence,info','o-')