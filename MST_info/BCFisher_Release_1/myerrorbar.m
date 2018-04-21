function myerrorbar(XX,YY,EE,col,shade)
%% Draw error bars
% Copyright (c) 2015, Ingmar Kanitscheider and Ruben Coen Cagli. 
% All rights reserved.
% See the file LICENSE for licensing information.

if(~exist('shade'))
    shade=0;
end

if(~shade)
    if (size(EE,1)==1 || size(EE,2)==1) % symmetric error
        XX=XX(:);
        EE=EE(:);
        YY=YY(:);
        for i=1:numel(XX)
            plot(XX([i i]),[YY(i)-EE(i) YY(i)+EE(i)],'-','Color',col);
        end
    else % non-symmetric c.i.
        for i=1:numel(XX)
            plot(XX([i i]),[EE(i,1) EE(i,2)],'-','Color',col);
        end
    end
else
    if (size(EE,1)==1 || size(EE,2)==1) % symmetric error
        XX=XX(:);
        EE=EE(:);
        YY=YY(:);
        h=fill(XX([1:end end:-1:1]),[YY(1:end)-EE(1:end); YY(end:-1:1)+EE(end:-1:1)]',col);
    else % non-symmetric c.i.
        h=fill(XX([1:end end:-1:1]),[EE(1:end,1); EE(end:-1:1,2)],col);
    end
    set(h,'edgecolor',col);
end

end