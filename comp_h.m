
clear
load w_ole_out

nu_w = length(w_ole_out)-1;

aux_w = w_ole_out;
clear w_ole_out;

%% Remove the last weight, the one corresponding to the bias
w_ole_out = aux_w(1:nu_w);


for j=1:nu_w
    h(j) = sum(w_ole_out(1:j));
end

middle = round(nu_w/2);

aux_h = h;
h(nu_w-middle+1:nu_w) = aux_h(1:middle);
h(1:nu_w-middle) = aux_h(middle+1:nu_w);

h(1:50)=h(1:50)+h(100:-1:51)/2;
h(51:100)=h(50:-1:1);
h=h-min(h);


h_mat(1,:) = h;
for j=2:nu_w
   h_mat(j,1:j-1) = h(nu_w-j+2:nu_w);
   h_mat(j,j:nu_w) = h(1:nu_w-j+1);
end




subplot(221)
plot(h)

subplot(222)
plot(exp(h))

subplot(223) 
plot(w_ole_out)

subplot(224)
surfl(h_mat)

%save h h h_mat