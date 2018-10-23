function SetFigure(size)
% Make figure more beautiful. @HH2013

if nargin == 0
    size = 15;
end

% === Colors ===
set(gcf,'color','w'); % Color
% set(findall(gcf,'-property','fontsize'),'fontsize',size); % Font size

% === Fonts ===
% set(findall(gcf,'type','axes'),'fontsize',size)
set(findall(gcf,'-property','fontsize'),'fontsize',size);
 
% === Axis ===
set(findall(gcf,'tickdir','i'),'tickdir','o'); % Tick direction
set(findall(gcf,'type','axes'),'ticklength',[0.02 0],'LineWidth',1.5,'color','none'); % Tick length
% set(findall(gcf,'type','axes','linewidth',0.5,'-not','tag','legend'),'linewidth',15); % Tick width

set(findall(gcf,'type','axes','-not','tag','legend'),'box','off'); % Box off

% === Plots ===
try
    set(findobj('type','errorbar'),'capsize',10)
end