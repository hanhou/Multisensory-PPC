function distance = dist(param,data)

amp= param(1);
K= param(2);
peak= param(3);
offset = param(4);

aux = size(data);
nu_pos = aux(2);
pos = [0:180/nu_pos:179.9999];

hill= amp*exp(K*(cos((pos-peak)/180*2*pi)-1))+offset;

distance = sum((data-hill).^2);

