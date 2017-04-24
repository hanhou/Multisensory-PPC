clear
load loglik32
load loglik64
load loglik128
load loglik256
load loglik512
load loglik128512
load loglik256512

allloglik = [loglik128512' loglik32' loglik64' loglik128' loglik256' loglik512' loglik256512'];
plot([0 50 100 150 200],allloglik,'-s');