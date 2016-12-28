function distance = dist2(param,in,out)

alpha= param(1);
amp= param(2);
curve = amp*log(1+exp(in/alpha));


distance = sum((out-curve).^2);

