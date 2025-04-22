function [colormax,colormin] = perclip(dat,perc)
%caxis([colormin colormax])
%dat : data
%perc : clip percent of the number of data
[m,n] = size(dat);
np = round(m*n*perc);
dat = sort(dat(:));
colormax = dat(m*n-np);
colormin = dat(np);

end