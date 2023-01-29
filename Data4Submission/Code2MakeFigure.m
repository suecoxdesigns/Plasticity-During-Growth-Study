

%% add work and av power into the structure


%% Jonas Plots 4 paper
formatOut = 'dd mmm yyyy';
clear

load('Vo2JumpData91417.mat')

%load(['MatlabCode\Data\ElasticData',datestr(now, formatOut),'.mat'])
load(['ElasticData25 Nov 2020.mat'])
d2 = es;
de = d2(([d2.group]=='E')|([d2.group]=='D'))

 for i=1:length(de)
     for j=1:length(de(i).PeNfl)
     de(i).PeNfl(j).avNfl = nanmean([de(i).PeNfl(j).MgNfl,de(i).PeNfl(j).MgNfl]);
   %  de(i).fMaxNFL = de(i).fMaxLg+de(i).fMaxMg + 
   
     end
     de(i).avOFL = nanmean([de(i).LgAvOfl, de(i).MgAvOfl]);
 end


 subplot(3,3,[1,2,4,5,7,8]); 


 
for i=1:length(de)
    % plot PE vs nFL
    de(i).fmax = de(i).fMaxLg+de(i).fMaxMg;
    thick = 10*((de(i).fmax-nanmin([de.fmax]))/(nanmax([de.fmax]-nanmin([de.fmax])))+1);
   if de(i).group=='D'
        col = [0/255,76/255,153/255];
    else
        col = [130/255,130/255,130/255];
    end

    if ~isempty(de(i).PeNfl)
        %col=(de(i).fMaxLg-nanmin([de.fMaxLg]))/(nanmax([de.fMaxLg])-nanmin([de.fMaxLg]));
        scatter1 =scatter([de(i).PeNfl.avNfl ],[de(i).PeNfl.PE],thick,'o','MarkerEdgeColor','k','MarkerFaceColor',col)
        scatter1.MarkerFaceAlpha = 0.2;
        scatter1.MarkerEdgeAlpha = 0.2;
        hold on
       % plot(de(i).MgLenA0c, de(i).PEc,'o','Color',col,'MarkerFaceColor','r')
    end
end

for i=1:length(de)
   plot(de(i).MgLenA0c, de(i).PEc,'o','Color',col,'MarkerFaceColor','r')  
end
xlabel('Av.Gastroc Passive N.Fiber Length')
ylabel('Elastic Strain Energy Stored (Max PE) J')

de= de(~cellfun(@isempty,{de.PE}))
set(gca, 'FontName', 'Times New Roman','FontSize',14)
hold off

ax1 = subplot(3,3,3)
for i=1:length(de)
   if de(i).group=='D'
        col = [67,90,138]/255;
    else
        col = [130,130,130]/255;
    end
    plot(de(i).fMaxLg+de(i).fMaxMg,de(i).PE,'o','Color',col,'MarkerFaceColor',col,'MarkerSize',5)
    hold on
   
end
xlabel('Max Iso Gastroc Force N')
ylabel('Maximum PE J')
ylim([0 0.55])

% fmax = [de.fMaxLg] + [de.fMaxMg];
% PE = [de.PE]
% P = polyfit(fmax,PE,1);
% %
% x0 = min(fmax) ; x1 = max(fmax) ;
% xi = linspace(x0,x1) ;
% yi = P(1)*xi+P(2);
% hold on
% plot(xi,yi,'k') ;


hold off
set(gca, 'FontName', 'Times New Roman','FontSize',14)

subplot(3,3,6)
for i=1:length(de)
    if de(i).group=='D'
        col = [67,90,138]/255;
    else
        col = [130,130,130]/255;
    end
    plot(de(i).tendonK/1000,de(i).PE,'o','Color',col,'MarkerFaceColor',col,'MarkerSize',5)
    hold on
end
xlabel('Tendon Stiffness kN/m')
ylabel('Maximum PE J')
xlim([20 80])
ylim([0 0.55])
 hold off
 
set(gca, 'FontName', 'Times New Roman','FontSize',14)


subplot(3,3,9)
for i=1:length(de)
    if de(i).group=='D'
        col = [67,90,138]/255;
    else
        col = [130,130,130]/255;
    end 
    if ~isempty(de(i).PE)
    plot(de(i).avOFL *100,de(i).PE,'o','Color',col,'MarkerFaceColor',col,'MarkerSize',5)
    hold on
   
    end
end
xlabel('Av. Muscle OFL ')
ylabel('Maximum PE J')
ylim([0 0.6])
xlim([20 36])
 hold off
 
set(gca, 'FontName', 'Times New Roman','FontSize',14)

% print('PaperFiles/PEvsNflDec2020','-dsvg','-r300');
% print('PaperFiles/Figure4Dec2020','-dpng','-r300');
hold off

