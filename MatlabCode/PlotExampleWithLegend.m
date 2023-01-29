x = 1:50;
y = 1:4:200;
y2 = y*2;
y3 = y*4;

plot(x,y,'s','Color','g')
fit1 = polyfit(x,y,2);
fity1 = polyval(fit1,x);
hold on
h1 = plot(x,fity1,'Color','g');
fitTxt = strcat('Copper: y=',num2str(fit1(2)),'x+',num2str(fit1(1)));

plot(x,y2,'s','Color','b')
fit2 = polyfit(x,y2,2);
fity2 = polyval(fit2,x);
hold on
h2 = plot(x,fity2,'Color','b');
fit2Txt = strcat('Carbon: y=',num2str(fit2(2)),'x+',num2str(fit2(1)));

plot(x,y3,'s','Color','k')
fit3 = polyfit(x,y3,2);
fity3 = polyval(fit3,x);
hold on
h3 = plot(x,fity3,'Color','k');
fit3Txt = strcat('Lead: y=',num2str(fit3(2)),'x+',num2str(fit3(1)));

legend([h1,h2,h3],{fitTxt,fit2Txt,fit3Txt},'Location','northwest')