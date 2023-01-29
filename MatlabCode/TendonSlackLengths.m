% this program will calculate the tendon slack length for the exercise and
% disuse birds at 6 months

% run the code to compile other data first
CombineVo2AndJumping2019_3
clear


formatOut = 'dd mmm yyyy';
load(['Data\MuscleTendonFunctionalData',datestr(now, formatOut),'.mat'])


% first read in muscle properties spreadsheet
%load('MuscleTendonFunctionalData22 Sep 2019.mat')
%load('C:\Users\smc288\Box Sync\Suzanne_Cox\SuzanneGF\Botox and Growth\SixMonthFunctionalPaper\Analyses\MuscleTendonFunctionalData09 Oct 2019.mat')

% link kavya's data to the 


musclePassiveForce = [0.8,0.998,0.999,1:0.1:1.6,1.601,1.602,2;0,0,0,0,0.035,0.12,0.26,0.55,1.17,2.0,2,2,2];
%individual muscle passive force scaled to fmax

% Step 1 - calculate muscle max force along tendon
%Fmax = (mass/(1.06*OFL*100))*30;
for i=1:length(d)
    if ~(isempty(d(i).LgAvOfl) || isempty(d(i).MgAvOfl) )
        
        % step 1: calculate Fmax for Lg and MG
        d(i).fMaxLg = (d(i).LG/(100*1.06*d(i).LgAvOfl/1000))*30;  % LG in g, Ofl in meters   % ask about pennation angle
        d(i).fMaxMg = (d(i).MG/(100*1.06*d(i).MgAvOfl/1000))*30;
        d(i).nLgFl = d(i).avLgFl/d(i).LgAvOfl;
        d(i).nMgFl = d(i).avMgFl/d(i).MgAvOfl;
         %density = 1060 kg/m^3   %specific tension of 3e5 N/M^2
        % CSA = Mmusc/(density*OFL)
        
        % Step 2 - calculate pennation angle at optimum
        %Optimal Pennation angle
        %oPA = asin((((sin(PennationAngle))*FiberLength))/OFL;
        
        
        
        d(i).optPenAngLg = asind(sind(d(i).lgPenAngle)*(d(i).avLgFl/d(i).LgAvOfl));
        avMgPenAngle = mean([d(i).mgPPenAngle,d(i).mgAPenAngle]);
        d(i).optPenAngMg = asind(sind(avMgPenAngle)*(d(i).avMgFl/d(i).MgAvOfl));
        d(i).fMaxLgT = d(i).fMaxLg*cosd(d(i).optPenAngLg);  % force along tendon
        d(i).fMaxMgT = d(i).fMaxMg*cosd(d(i).optPenAngMg);
        
        %Step 3 calculate fMax along tendon
        % fmaxT = fmax*cos(optPenAng)
        if (isempty(d(i).optPenAngLg) || isempty(d(i).optPenAngMg))
            continue
        end

        % step 4 for each muscle - scale the passive f-l curve
        
        
        d(i).pasMuscleLg = [musclePassiveForce(1,:);musclePassiveForce(2,:)*d(i).fMaxLgT];
        d(i).pasMuscleMg = [musclePassiveForce(1,:);musclePassiveForce(2,:)*d(i).fMaxMgT];
        
        
        % Step 5 Calculate the passive muscle force at the normalized fiber length
        
        
        d(i).nLgFl = d(i).avLgFl/d(i).LgAvOfl;
        d(i).nMgFl = d(i).avMgFl/d(i).MgAvOfl;
        if ~isnan(d(i).fMaxLgT)
            d(i).pasForceLg = spline(d(i).pasMuscleLg(1,:),d(i).pasMuscleLg(2,:),d(i).nLgFl);
            if d(i).pasForceLg<0
                d(i).pasForceLg = 0;
            end
        end
        if ~(isnan(d(i).fMaxMgT)|| isempty(d(i).tendonCrossAv))
            d(i).pasForceMg = spline(d(i).pasMuscleLg(1,:),d(i).pasMuscleLg(2,:),d(i).nMgFl);
            if d(i).pasForceMg<0
                d(i).pasForceMg = 0;
            end
            d(i).totPassForceLgMg = d(i).pasForceMg + d(i).pasForceLg;
            d(i).tendPassStress = d(i).totPassForceLgMg/(d(i).tendonCrossAv/(1000^2));
        end
    end
end

%% load the tendon data

%Trim and resample the data
%TendonShifting
 
      clearvars -except d  
    % First load Kavya's data
    load('Data\Tendon4.mat')
    % combine the data for each animal and fit a spline
    [d.AniNum] = deal(d.aniNum);
    d = rmfield(d,'aniNum');
    Tendon = ShareData;
    clear ShareData
    d([d.AniNum]==797).AniNum = 796;
    Tendon(11).ShiftedData.Cycle(11).DispShifted =  Tendon(11).ShiftedData.Cycle(11).DispShifted(1:end-1);
    Tendon(11).ShiftedData.Cycle(11).ForceShifted =  Tendon(11).ShiftedData.Cycle(11).ForceShifted(1:end-1);
    
    for i=1:length(Tendon)
        disp = [];
        force = [];
        
        dIdx = find([d.AniNum] == Tendon(i).BirdNumber);
        if isempty(dIdx)
            dIdx = length(d)+1;
            d(dIdx).AniNum = Tendon(i).BirdNumber;
            continue
        end
        if d(dIdx).AniNum == 761
            d(dIdx).AniNum
        end
        [maxVal,idxMax] = max( mean([Tendon(i).ShiftedData.Cycle.DispShifted],2));
        meanDist =  mean([Tendon(i).ShiftedData.Cycle.DispShifted],2);
        meanForce = mean([Tendon(i).ShiftedData.Cycle.ForceShifted],2);
        Tendon(i).cDist = meanDist(1:round(idxMax*0.9));
        Tendon(i).cForce =meanForce(1:round(idxMax*0.9));
        Tendon(i).cStrain = [Tendon(i).cDist]./Tendon(i).TendonLength;
        
        if ~(isempty(Tendon(i).cStrain))
            clear disp force stress modulous newStrain newF
            strain = Tendon(i).cStrain(1:round(idxMax*0.9));
            forceM = Tendon(i).cForce(1:round(idxMax*0.9));
   
            if (~isempty(d(dIdx).tendonCrossMin)&&~isnan(d(dIdx).tendonCrossMin))
                stress = Tendon(i).cForce/(d(dIdx).tendonCrossAv/(1000^2));
                  mod = polyfit(strain(end-60:end),stress(end-60:end),1);
                modulous =mod(1);
                d(dIdx).tendonE = modulous;
            end
            tendonPassiveStrain = spline(forceM,strain,d(dIdx).totPassForceLgMg);
            d(dIdx).tendonPassiveDispMm = tendonPassiveStrain*((Tendon(i).TendonLength));
            d(dIdx).tendonPassiveStrain =  tendonPassiveStrain;
            d(dIdx).tendonSampleLength = Tendon(i).TendonLength;
            p = polyfit(Tendon(i).cDist(end-50:end)/1000,forceM(end-50:end),1)
            k = p(1);
            d(dIdx).tendonK = k;
            
            newStrain = strain;
            newF=Tendon(i).cForce;
            
            %find the slope of the last two points - extrapolate out to
            %high strain
            slope = (newF(end)-newF(end-50))/(newStrain(end)-newStrain(end-50));
            b = newF(end)-slope*strain(end);
            newYs = slope.*([9.2 9.201 9.202 12])+b;
            tendonForceN =[newF',newYs];
            
            
            % plot(x,y)
            %[curve, goodness, output] = fit(x,y,'smoothingspline','SmoothingParam',0.8);
            %    plot(xs,ys)
            spline1 = spap2(1,4,newStrain, newF);
            newx = [newStrain(1):0.002:max(newStrain)];
            newy = fnval(spline1,newx);
            
            
            
            d(dIdx).tendonForce = [newy,newYs];
            d(dIdx).tendonStrain= [newx,9.2 9.201 9.202 12];
            
            %             plot(newStrain,newF)
            %             hold all
        end
        
        %, 9.2 9.201 9.202 12] tendon strain extra points
    end
    %%
  
  % Last step - calculate the tendon slack length
  
  % MTU length = MuscleLength + TendonSlackLength +
  % TendonSlackLength*TendonStrain
  
  for i=1:length(d)
      if ~isnan(d(i).tendonPassiveDispMm)
          avMgPenAngle = mean([d(i).mgPPenAngle,d(i).mgAPenAngle]);
          LgFiberLengthAlongTendon = d(i).avLgFl*cosd(d(i).lgPenAngle);
          MgFiberLengthAlongTendon = d(i).avMgFl*cosd(avMgPenAngle);
          
          d(i).tendonSlackLengthLGmm = (d(i).mtuLengLG - LgFiberLengthAlongTendon)/(1+d(dIdx).tendonPassiveStrain);
          d(i).tendonSlackLengthMGmm = (d(i).mtuLengMG - MgFiberLengthAlongTendon)/(1+d(dIdx).tendonPassiveStrain);
          
      end
  end
  
  save('Data\tendonSlackLengthsJan2020.mat','d')

%% Add in moment arm data

load('Data\ZanneMaDataV5.mat')

% re-order and arrange the moment arm data
MA = struct('angles',{},'ma',{},'aniNum',{})
birds = unique([mas.BirdID]);
cats = categories([mas.BirdID]);


for i=1:length(birds)
    MA(i).aniNum = str2num(cats{i}(end-2:end));
    if MA(i).aniNum == 795
        i
    end
    sub = mas([mas.BirdID] == birds(i));
    ang = unique([sub.Angle]);
    count=1;
    for j=1:length(ang)
        if ~isnan(ang(j))
            sub2 = sub([sub.Angle]==ang(j));
            tempMA = nanmean([sub2.MA]);
            if ~isnan(tempMA)
                MA(i).ma(count) =tempMA ;
                MA(i).angles(count) = ang(j);
                count=count+1;
            end
        end
    end
    %   plot([sub.Angle],[sub.MA],'*');
%   hold on
%   plot([MA(i).angles],[MA(i).ma]);
%   pause
end
  

for i=1:length(MA)
    dIdx = find([d.AniNum] == MA(i).aniNum);
    if isempty(dIdx)
        if ~isnan(MA(i).aniNum)
        dIdx = length(d)+1;
        d(dIdx).AniNum = MA(i).aniNum;
        else
            continue
        end
    end
    d(dIdx).avMa = nanmean(MA(i).ma)/10;
    d(dIdx).ma = MA(i).ma/1000;
    d(dIdx).maAngles = MA(i).angles;
    clear dIdx
    
end

%save new MA data
save('Data\MomentArmData2.mat','MA')


%% Add in leg length data
load('Data\LegLengths.mat')
load('Data\genericLimbLengths.mat')

birds = unique([ll.BirdNum]);
for i=1:length(birds)
    subset = ll([ll.BirdNum]==birds(i));
    if birds(i)==711
        birds(i)=796;
    end
    if birds(i)==781
        birds(i)=785;
    end
    if birds(i) == 790
        birds(i)=760;
    end
    
    femur = nanmean([subset.femur]);
    tibia = nanmean([subset.tibia]);
    tmt = nanmean([subset.TMT]);
    toe = nanmean([subset.Toe]);
    
    idx = find([d.AniNum] == birds(i));
    if ~isempty(idx)
        d(idx).femur = femur;
        d(idx).femurN = femur/lgts.femur;
        d(idx).tibia = tibia;
        d(idx).tibiaN = tibia/lgts.tibia;
        d(idx).tmt = tmt;
        d(idx).tmtN = tmt/lgts.tmp;
        d(idx).toe = toe;
        d(idx).toeN = toe/lgts.toes;
    else
       idx = length(d)+1;
        d(idx).femur = femur;
        d(idx).femurN = femur/lgts.femur;
        d(idx).tibia = tibia;
        d(idx).tibiaN = tibia/lgts.tibia;
        d(idx).tmt = tmt;
        d(idx).tmtN = tmt/lgts.tmt;
        d(idx).toe = toe;
        d(idx).toeN = toe/lgts.toe;
        d(idx).AniNum = birds(i);
    end
    clear femur idx tibia tmt toe femurN tibiaN toeN tmtN
end

%% convert to units used in OSimm

for i=1:length(d)
    d(i).optPenAngMg = deg2rad(d(i).optPenAngMg);
    d(i).optPenAngLg = deg2rad(d(i).optPenAngLg);
    d(i).MgAvOfl = (d(i).MgAvOfl)/100;
    d(i).LgAvOfl = (d(i).LgAvOfl)/100;
    d(i).maxFAlongTendonMG = d(i).fMaxMg*cosd(d(i).optPenAngMg);
    d(i).maxFAlongTendonLG = d(i).fMaxLg*cosd(d(i).optPenAngLg);
    if ~isempty(d(i).fMaxLg)
        d(i).tendonForceNLG = d(i).tendonForce/d(i).maxFAlongTendonLG;
    end
    if ~isempty(d(i).fMaxMg)
        d(i).tendonForceNMG = d(i).tendonForce/d(i).maxFAlongTendonMG;
    end
end
d = orderfields(d);
d([d.AniNum]==796).avMa = d([d.AniNum]==797).avMa;
d([d.AniNum]==796).maAngles = d([d.AniNum]==797).maAngles;
d([d.AniNum]==796).ma = d([d.AniNum]==797).ma;

%% Simplify data to put into Osim and R
d([d.AniNum]==758).group = 'B';
d= orderfields(d);
dn=fieldnames(d);
%save tendon min CSA, tendon modulous,MG, LG,
%'BmdGcm2','PerFat','tendonCrossMin'

tokeep = {'AniNum','group','mass','mtuLengLG','avLgFl','avMgFl','nLgFl','nMgFl','avRvO2','mtuLengMG','tendonSampleLength','sumExtensors','maxJump','ma','maxPower','maxV','maAngles','MGLength','LGLength','tendonK','maxFAlongTendonLG','tendonForceNLG','tendonForceNMG','maxFAlongTendonMG','femurN','femur','tibia','tibiaN','tmt','tmtN','toe','toeN','fMaxLgT','fMaxMgT','tendonCrossMin','tendonCrossAv','tendonE','PerFat','BmdGcm2','MgAvOfl','LG','MG','avMa','LgAvOfl','optPenAngLg','optPenAngMg','fMaxLg','fMaxMg','tendonForce','tendonStrain','tendonSlackLengthMGmm','tendonSlackLengthLGmm'};
toRemove = dn(~ismember(dn,tokeep));

d2 = rmfield(d,[toRemove]);



count = length(d2);
while count>0
    if (isempty(d2(count).tendonSlackLengthMGmm)&& isempty(d2(count).avMa))
        d2(count) = [];
        if count>length(d2)
            count=count-1;
        end
    else
        count=count-1;
    end
end
formatOut = 'dd mmm yyyy';
save(['Data\MorphData',datestr(now, formatOut),'.mat'],'d2')
d3 = rmfield(d2,{'tendonStrain','tendonForce','tendonForceNLG','tendonForceNMG','ma','maAngles'});
Rdata2 = struct2table(d3);
writetable(Rdata2,['Data\MorphData',datestr(now, formatOut),'.csv']);

clearvars -except d2