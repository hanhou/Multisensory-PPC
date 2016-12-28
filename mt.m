clear

r_spont_MT = 20;
b_pref = 0.4;
b_null = -0.2;
K_in = 4;
N_in = 100;
pos_in=[0:180/N_in:179.9999];
delta_t = 1e-3;

coherence =0;
for j=1:6
   max_rate_in =r_spont_MT + b_pref*coherence;
   b_in = r_spont_MT + b_null*coherence;
   MT(j,:)=((max_rate_in-b_in)*exp(K_in*(cos((pos_in'-90)/180*2*pi)-1))+b_in); 
   coherence= 3*2^(j-1);
end

plot([0:100/N_in:99.9999],MT');
xlabel('Stimulus');
ylabel('Firing rate (Hz)');