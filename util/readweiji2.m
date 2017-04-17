clear

load exp_data512

[junk nu_trial] = size(spike_time50);

spikes2_out50 = zeros(nu_trial,401); 
spikes2_out100 = zeros(nu_trial,401); 

for k=1:nu_trial
   nu_spikes = length(spike_time50{k}); 
   for l=1:nu_spikes
       spikes2_out50(k,spike_time50{k}(l))=1;
   end   
   nu_spikes = length(spike_time100{k}); 
   for l=1:nu_spikes
       spikes2_out100(k,spike_time100{k}(l))=1;
   end
end



subplot(221)
plot(sum(spikes2_out50),'r')
hold on
plot(sum(spikes2_out100))
hold off

subplot(222)
plot(sum(proba_out50),'r')
hold on
plot(sum(proba_out100))
hold off


subplot(223)
plot(sum(spikes2_out50),'r')
hold on
plot(sum(proba_out50)/10^3,'r')
plot(sum(proba_out100)/10^3)
plot(sum(spikes2_out100))
hold off
axis([0 300 0 120]);


subplot(224)
plot(sum(spikes_out100))
hold on
plot(sum(spikes_out50),'r')
hold off
