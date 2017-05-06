function [h_column_wise, h_row_wise] = tight_subplot(Nh, Nw, gap, marg_h, marg_w, ratio_h, ratio_w)

% tight_subplot creates "subplot" axes with adjustable gaps and margins
%
% ha = tight_subplot(Nh, Nw, gap, marg_h, marg_w, ratio_h, ratio_w)
%
%   in:  Nh      number of axes in hight (vertical direction)
%        Nw      number of axes in width (horizontaldirection)
%        gap     gaps between the axes in normalized units (0...1)
%                   or [gap_h gap_w] for different gaps in height and width 
%        marg_h  margins in height in normalized units (0...1)
%                   or [upper lower] for different lower and upper margins 
%        marg_w  margins in width in normalized units (0...1)
%                   or [left right] for different left and right margins 
%
%  out:  [h_column_wise, h_row_wise]     array of handles of the axes objects
%                   starting from upper left corner, going column-wise and row-wise, respectively
%
%  Example: ha = tight_subplot(3,2,[.01 .03],[.1 .01],[.01 .01])
%           for ii = 1:6; axes(ha(ii)); plot(randn(10,ii)); end
%           set(ha(1:4),'XTickLabel',''); set(ha,'YTickLabel','')

% Pekka Kumpulainen 20.6.2010   @tut.fi
% Tampere University of Technology / Automation Science and Engineering

% Improved by HH20140526


if nargin<3 || isempty(gap); gap = .02; end
if nargin<4 || isempty(marg_h); marg_h = .05; end
if nargin<5 || isempty(marg_w); marg_w = .05; end
if nargin<6 || isempty(ratio_h); ratio_h = ones(1,Nh); end
if nargin<7 || isempty(ratio_w); ratio_w = ones(1,Nw); end

if numel(gap)==1; 
    gap = [gap gap];
end
if numel(marg_w)==1; 
    marg_w = [marg_w marg_w];
end
if numel(marg_h)==1; 
    marg_h = [marg_h marg_h];
end

% HH20140526
% axh = (1-sum(marg_h)-(Nh-1)*gap(1))/Nh; 
% axw = (1-sum(marg_w)-(Nw-1)*gap(2))/Nw;

axh = (1-sum(marg_h)-(Nh-1)*gap(1)) * ratio_h/ sum(ratio_h); 
axw = (1-sum(marg_w)-(Nw-1)*gap(2)) * ratio_w/ sum(ratio_w);

% py = 1-marg_h(2)-axh; 

h_column_wise = zeros(Nh*Nw,1);
h_row_wise = zeros(Nh*Nw,1);
% ii = 0;

% HH20140526
for ii = 1 : Nh*Nw
    [i,j] = ind2sub([Nh Nw],ii);
    
    px = marg_w(1) + gap(2) * (j-1) + sum(axw(1:j-1));
    
    py = marg_h(2) + gap(1) * (Nh - i) + sum(axh(i+1:end));
    h_column_wise(ii) = axes('Units','normalized', ...
        'Position',[px py axw(j) axh(i)]);
    
    h_row_wise(sub2ind([Nw,Nh],j,i)) = h_column_wise(ii);
    
        %'XTickLabel','', ...
       % 'YTickLabel',''
end




% 
% for ih = 1:Nh
%     px = marg_w(1);
%     
%     for ix = 1:Nw
%         ii = ii+1;
%         ha(ii) = axes('Units','normalized', ...
%             'Position',[px py axw() axh()], ...
%             'XTickLabel','', ...
%             'YTickLabel','');
%         px = px+axw+gap(2);
%     end
%     py = py-axh-gap(1);
% end

