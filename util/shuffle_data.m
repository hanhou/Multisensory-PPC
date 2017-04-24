      for j=1:N_out

         for m=1:2
         for n=1:nu_window
            spikes_count_shuff{m,n}(j,1:nu_trial) = ...
                spikes_count{m,n}(j,mod(j-1:nu_trial-1+(j-1),nu_trial)+1);    
            spikes_count_shuff{m,n}(j,nu_trial+1:nu_trial*2) = ...
                spikes_count{m,n}(j,mod(j-1:nu_trial-1+(j-1),nu_trial)+nu_trial+1);    
         end
         end
         
         spikes_out_shuff_trial11(j,1:nu_trial) = spikes_out_trial11(j,mod(j-1:nu_trial-1+(j-1),nu_trial)+1);    
         spikes_out_shuff_trial12(j,1:nu_trial) = spikes_out_trial12(j,mod(j-1:nu_trial-1+(j-1),nu_trial)+1);  
         spikes_out_shuff_trial11(j,nu_trial+1:nu_trial*2) = spikes_out_trial11(j,mod(j-1:nu_trial-1+(j-1),nu_trial)+nu_trial+1);    
         spikes_out_shuff_trial12(j,nu_trial+1:nu_trial*2) = spikes_out_trial12(j,mod(j-1:nu_trial-1+(j-1),nu_trial)+nu_trial+1); 
         
         spikes_out_shuff_trial21(j,1:nu_trial) = spikes_out_trial21(j,mod(j-1:nu_trial-1+(j-1),nu_trial)+1);    
         spikes_out_shuff_trial22(j,1:nu_trial) = spikes_out_trial22(j,mod(j-1:nu_trial-1+(j-1),nu_trial)+1);  
         spikes_out_shuff_trial21(j,nu_trial+1:nu_trial*2) = spikes_out_trial21(j,mod(j-1:nu_trial-1+(j-1),nu_trial)+nu_trial+1);    
         spikes_out_shuff_trial22(j,nu_trial+1:nu_trial*2) = spikes_out_trial22(j,mod(j-1:nu_trial-1+(j-1),nu_trial)+nu_trial+1);  

         spikes_out_shuff_trial31(j,1:nu_trial) = spikes_out_trial31(j,mod(j-1:nu_trial-1+(j-1),nu_trial)+1);    
         spikes_out_shuff_trial32(j,1:nu_trial) = spikes_out_trial32(j,mod(j-1:nu_trial-1+(j-1),nu_trial)+1);  
         spikes_out_shuff_trial31(j,nu_trial+1:nu_trial*2) = spikes_out_trial31(j,mod(j-1:nu_trial-1+(j-1),nu_trial)+nu_trial+1);    
         spikes_out_shuff_trial32(j,nu_trial+1:nu_trial*2) = spikes_out_trial32(j,mod(j-1:nu_trial-1+(j-1),nu_trial)+nu_trial+1);  
       
         spikes_out_shuff_trial41(j,1:nu_trial) = spikes_out_trial41(j,mod(j-1:nu_trial-1+(j-1),nu_trial)+1);    
         spikes_out_shuff_trial42(j,1:nu_trial) = spikes_out_trial42(j,mod(j-1:nu_trial-1+(j-1),nu_trial)+1);  
         spikes_out_shuff_trial41(j,nu_trial+1:nu_trial*2) = spikes_out_trial41(j,mod(j-1:nu_trial-1+(j-1),nu_trial)+nu_trial+1);    
         spikes_out_shuff_trial42(j,nu_trial+1:nu_trial*2) = spikes_out_trial42(j,mod(j-1:nu_trial-1+(j-1),nu_trial)+nu_trial+1); 
      end