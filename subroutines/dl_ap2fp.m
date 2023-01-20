function [fp] = yc_ap2fp(ap)
% yc_ap2fp: axis point to figure point (e.g., for consistent and easy annotation)
%
% INPUT
% ap: input point in axis (x,y)
%
% OUTPUT
% fp: output point in figure (x,y)
% 
% DEMO
% test/test_plot_yc_ap2fp.m
%
% TODO
% If you have time, make a summary of Matlab annotation strategies 
% (arrow, textarrow, ellipse, recntangles, etc)

apos=get(gca,'Position');%[x,y,wx,wy], wx: width in x
xm=xlim;    %xm=[xmin,xmax]
ym=ylim;    %ym=[ymin,ymax]

if string(get(gca,'XDir'))=='reverse' %from large x to small
fp(1)=(xm(2)-ap(1))/(xm(2)-xm(1)) * apos(3)+apos(1);%x coordinate in figure
else                          %from small x to large
fp(1)=(ap(1)-xm(1))/(xm(2)-xm(1)) * apos(3)+apos(1);%x coordinate in figure
end


if string(get(gca,'YDir'))=='reverse' %from large y to small
fp(2)=(ym(2)-ap(2))/(ym(2)-ym(1)) * apos(4)+apos(2);%y coordinate in figure
else                          %from small y to large
fp(2)=(ap(2)-ym(1))/(ym(2)-ym(1)) * apos(4)+apos(2);%y coordinate in figure 
end




return