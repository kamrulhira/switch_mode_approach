function [PEARSON] =Plot2figure(Xp,Yp,str,ERR3,sub,save_on)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
PEARSON=corr(Xp',Yp');
hfig=figure;
plot(Xp,Yp);
    h3=title(['Pearson Correlation: ',num2str(PEARSON)]);
    set(h3,'fontsize',20,'FontName','Times New Roman','fontweight','bold');
    opt.XLabel = 'Ground-Truth  of  Heart  Rate (BPM)';opt.YLabel = 'Estimates  of  Heart  Rate (BPM)'; %ylabel
    opt.XLim = [60,200 ];opt.YLim = [60,200 ];
opt.Colors = [ % three colors for three data set
    1,      0,       0;
    ];
% opt.LineWidth = [3, 3]; % three line widths
opt.LineStyle = {'None'}; % three line styles
opt.Markers = { 'o'};
% opt.Legend = {[ 'Pearson Correlation: ' num2str(PEARSON)]}; % legends
opt.BoxDim = [5,5]; %[width, height]

setPlotProp(opt,hfig);
if save_on==1
kk=[ 'Pearson.png'];
saveas(hfig,[pwd str '/' kk]);
kk=[ 'Pearson.eps'];
saveas(hfig,[pwd str '/' kk]);
kk=[ 'Pearson.tif'];
saveas(hfig,[pwd str '/' kk]);
kk=[ 'Pearson.jpg'];
saveas(hfig,[pwd str '/' kk]);
end
hfig=figure;plot((Xp+Yp)/2,Xp-Yp);hold on;
u=mean(Xp-Yp); 
SSD=mean(ERR3(sub)); 
LOA=[u+SSD*1.96 u-SSD*1.96];
zzz=u+ones(1,141)*SSD*1.96;
zzz2=u-ones(1,141)*SSD*1.96;
zxx=ones(1,141)*u;
xax=60:200;
plot(xax,zzz),plot(xax,zzz2),plot(xax,zxx)
opt.XLabel=('Average  of  Ground-Truth  and  Estimate  (BPM)');
opt.YLabel='Difference ';
opt.XLim = [60,200 ];opt.YLim = [min(Xp-Yp),20];
opt.Colors = [ % three colors for three data set
    1,      0,       0;
    0,      0,       0;
    0,      0,       0;
    0,      0,       0;
    ];
opt.LineStyle = {'None',':',':',''}; % three line styles
opt.Markers = { '+','','',''};
opt.BoxDim = [5,5]; %[width, height]
setPlotProp(opt,hfig);
FS=20;
h4=text(180,LOA(1)+1.5,'\mu+1.96\sigma');
set(h4,'fontsize',FS,'FontName','Times New Roman','fontweight','bold','color','blue');
h4=text(185,u+1.5,'\mu');
set(h4,'fontsize',FS,'FontName','Times New Roman','fontweight','bold','color','blue');
h4=text(180,LOA(2)+1.5,'\mu-1.96\sigma');
set(h4,'fontsize',FS,'FontName','Times New Roman','fontweight','bold','color','blue');
% set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0, 1, 1]);
if save_on==1
kk=[ 'Bland.png'];
saveas(hfig,[pwd str '/' kk]);
kk=[ 'Bland.eps'];
saveas(hfig,[pwd str '/' kk]);
kk=[ 'Bland.tif'];
saveas(hfig,[pwd str '/' kk]);
kk=[ 'Bland.jpg'];
saveas(hfig,[pwd str '/' kk]);
end
finnal=1;
end

