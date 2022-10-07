cWvl = 632.8;
nWvl = 1000;
wvl = linspace(cWvl-10,cWvl+10,nWvl);
resp = ones(1,nWvl);
ind = wvl>cWvl+5;resp(ind) = 0;
ind = wvl<cWvl-5;resp(ind) = 0;

figure; 
hax = axes;
hold on
for i = -5:5
    SP=cWvl+i;
    
    l=line([SP SP],get(hax,'YLim'),'LineStyle','--','LineWidth',3,'Color',[0 0 1]);
    
end
p1=plot(wvl,resp,'k-','LineWidth',3);grid('on');
xlabel('Wavelength [nm]','FontSize',24)
ylabel('Response','FontSize',24)
axis([cWvl-10 cWvl+10 -0.1 1.2])

legend([p1,l],{'Continuous Filter','Discrete Sample'},'FontSize',24)