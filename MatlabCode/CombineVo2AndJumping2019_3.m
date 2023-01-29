%%  Combine Vo2 data with jump data and fix mass info
% before you run this - run 'CompileFiberLengthData_afterRedo.m' in the old code
% folder under analyses

CompileFiberLengthData_afterRedo2

clear  
load('Data\Vo2DataOct2019b.mat')
load('Data\SixMonthMassData.mat')  %use these masses
load('Data\AllCompiledJumpData-March2019.mat')
load('Data\Vo2JumpData318 Oct 2019.mat')
for i=1:length(d)
    if isempty(d(i).jump)
        d(i).jumpMaxNorm = 1;
        % d(i).forceTime = sum(d(i).jumpMaxNorm)/d(i).trialTimeSec;
       % d(i).zImpulses = 0;
       % d(i).xyImpulses = [0,0,0];
        %d(i).ImpulseRatio = [0,0,0];
    else
        d(i).jumpMaxZero = [d(i).jump.maxJump];
        d(i).jumpMaxNorm = [d(i).jump.maxJump]./d(i).weight;
        d(i).zImpulses = [d(i).jump.impulse]./d(i).weight;

        for j=1:length(d(i).jump)
            d(i).jump(j).jumpMaxNorm = d(i).jump(j).maxJump./d(i).weight;
            d(i).jump(j).impulse2 = d(i).jump(j).impulseStart2StandZ+d(i).jump(j).impulseStand2endZ;
        end
       d(i).zImpulses2 = [d(i).jump.impulse2]./d(i).weight;
       d(i).nJumps =  length(d(i).jump);
    end
end




d = rmfield(d,{'MaxJump','JumpMaxZ','JumpTimes','JumpMaxZZeroed','JumpMaxZNorm','JumpMaxZnum'});
formatOut = 'dd mmm yyyy';
save(['Data\CompiledJumpDataA,',datestr(now, formatOut),'.mat'],'d')
%boxplot([d.MaxJump],{d.Treat})


%% Combining data from multiple trials for each bird
clear
formatOut = 'dd mmm yyyy';
load(['Data\CompiledJumpDataA,',datestr(now, formatOut),'.mat'])
d2 = d;
for i=1:length(d)
    count = 1;
    for j=1:length(d(i).jump)
        if length(d(i).jump(j).jumpTLong)>=length(d(i).jump(j).jumpT)
            d2(i).jump(count) = d(i).jump(j);
            count=count+1;
        end
    end
end
d = d2;
aniNums = unique({d.AniNumS});
d2 = struct('AniNum',{},'Treat',{},'dates',{},'weight',{},'AniNumS',{});
count=1;
for i=1:length(aniNums)
    x1 = aniNums(i);
    y={d.AniNumS};
    substrfind = @(x,y) ~cellfun(@isempty,strfind(y,x));
    logicalArray1 = substrfind(x1,y);
    clear sublist filedates
    sublist= d(logicalArray1);
    if strcmp(sublist(1).Treat,'A')
        %before group
        subsetb = sublist([sublist.group]=='b');
        d2=extractJump2019c(subsetb,count,d2);
        count = count+1;
        %after group - subset further by date to get effect of time
        subseta = sublist([sublist.group]=='a');
        %divide up into different dates
        dates = unique({subseta.date2});
        
        for k=1:length(dates)
            s = subseta(strcmp({subseta.date2},dates(k)));
            [s.Treat]= deal(['A',num2str(k)]);
            d2=extractJump2019c(s,count,d2);
            count = count+1;
        end
    else
        %redo calculations with average mass
        d2=extractJump2019c(sublist,count,d2);
        count = count+1;
    end
end


% Subset data to Adult vs Growth studies
x1 = 'A';
y={d2.Treat};
substrfind = @(x,y) ~cellfun(@isempty,strfind(y,x));
logicalArray1 = substrfind(x1,y);
adult = d2(logicalArray1);
growth = d2(~logicalArray1);

x1 = 'B';
y={growth.Treat};
substrfind = @(x,y) ~cellfun(@isempty,strfind(y,x));
logicalArray1 = substrfind(x1,y);
d6 = growth(~logicalArray1);


%% combine jumping and VO2 data
load('Data\Vo2Data10-14-18ByBird.mat')
load('Data\SixMonthMassData.mat') 
sNames = unique([s.aniNum]);
dNames = unique([growth.AniNum]);
AniNums = unique([sNames,dNames]);
fNames = fieldnames(s);
d = struct('aniNum',{});



for j = 1:length(AniNums)
    d(j).aniNum = AniNums(j);
    for i=1:length(s)
        if (s(i).aniNum == AniNums(j))
            for f = 1:length(fNames)
            d(j).(fNames{f}) = s(i).([fNames{f}]);
            end
        end
    end
    
    for i=1:length(m)
        if (d(j).aniNum == m(i).AniNum)
            d(j).mass = m(i).Mass;
        end
    end
    for i=1:length(dNames)
        if growth(i).AniNum == AniNums(j)
            d(j).group = growth(i).Treat;
            d(j).meanJumpTimeTop3 = growth(i).meanJumpTimeTop3;
            d(j).maxJump = growth(i).MaxJumpTot;
            d(j).maxImpulseNonNorm = growth(i).maxImpulseNonZeroed;
            d(j).avTop3Jumps = growth(i).Avoftop3MaxJumps;
            d(j).maxImpulse = growth(i).maxImpulse;       
            d(j).jumpPsec  = growth(i).jumpPsec;
             d(j).nJumps = growth(i).nJumps;
             d(j).maxV = growth(i).maxV;
             d(j).CoefVarAv = growth(i).CoefVarAv;
             d(j).AvTop3Impulses = growth(i).AvTop3Impulses;
             d(j).jumpTimeTopImpulse = growth(i).jumpTimeTopImpulse;
            d(j).maxPower = growth(i).maxPower;
             
        end
        if (isfield(d,'mass') && ~isempty(d(j).avRvO2))
            d(j).navRVO2 = d(j).avRvO2/d(j).mass;
            d(j).navSVO2 = d(j).avSvO2/d(j).mass;
        end
    end
    if d(j).aniNum == 792
        d(j).group = 'D';
    elseif d(j).aniNum == 781
        d(j).group = 'E';
    end
end




formatOut = 'dd mmm yyyy';
save(['Data\VO2DataNormalized',datestr(now, formatOut),'.mat'],'d');
%% combine with muscle mass data
clear
% Setup the Import Options
opts = spreadsheetImportOptions("NumVariables", 35);

% Specify sheet and range
opts.Sheet = "DissectionMuscleMasses6MonthC";
opts.DataRange = "A2:AI26";

% Specify column names and types
opts.VariableNames = ["aniNum", "R_legmass", "L_legmass", "ILPO", "MG", "LG", "ILPR", "IC", "FL", "IG", "SupDF2", "SupDF3", "DF4", "TC", "DFpostXkneeplant", "FHL", "FDLfb", "FTD", "EDL", "IF", "Ambiens", "FTM", "FTICombineWithFTLP", "IIiof", "ITC", "ITCR", "ITM", "FCLP", "FCLA", "CF", "PIFM", "PIFL", "FCM", "ISHf", "OM"];
opts.VariableTypes = ["double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double"];

% Import the data
DissectionMuscleMasses6MonthD = readtable("C:\Users\smc288\Box\Suzanne_Cox\SuzanneGF\Botox and Growth\SixMonthFunctionalPaper\Analyses\Data\DissectionMuscleMasses6MonthD.xlsx", opts, "UseExcel", false);

% Clear temporary variables
dis = table2struct(DissectionMuscleMasses6MonthD);
save('C:\Users\smc288\Box\Suzanne_Cox\SuzanneGF\Botox and Growth\SixMonthFunctionalPaper\Analyses\Data\DissectionMuscleMasses6Month.mat','dis')
clear opts

%%
clear
load ('C:\Users\smc288\Box\Suzanne_Cox\SuzanneGF\Botox and Growth\SixMonthFunctionalPaper\Analyses\Data\DissectionMuscleMasses6Month.mat')
%load('VO2Data4-24-19Normalized.mat');
formatOut = 'dd mmm yyyy';
load(['Data\VO2DataNormalized',datestr(now, formatOut),'.mat']);

mFields = fieldnames(dis);
mFields = mFields(2:end);
flexors = {'IF','ITC','TC','ITM','EDL','ITCR','R_legmass','L_legmass'};
Gs = {'MG','LG'};
sumExtensor = 0;
for i=1:length(d)
    for j = 1:length(dis)
        if dis(j).aniNum==d(i).aniNum
           % d(i).sumExtensors = 0;
            d(i).rLegMass = dis(j).R_legmass;
            d(i).lLegMass = dis(j).L_legmass;
            for m = 2:length(mFields)
                d(i).(mFields{m}) = dis(j).(mFields{m});              
            end
        end  
    end
end


%break up into groups - botox - exercise - disuse
% fill in NAN's in muscle data with mean of group value for that muscle
subB = d([d.group]=='B');
subD = d([d.group]=='D');
subE = d([d.group]=='E');

d(1).extensors = 0;
d(1).flexors = 0;
d(1).Gs = 0;
for i=1:length(d)
    countE = 0;
    countF = 0;
    countGs  =0;
    for j=3:length(mFields)
        if ~isempty (d(i).rLegMass) % only do this for birds that we have muscle masses for
            if (~isnumeric(d(i).(mFields{j})))% if it is empty or nan - fill in with average
                switch d(i).group
                    case 'B'
                        d(i).(mFields{j}) = nanmean([subB.(mFields{j})]);
                    case 'D'
                        d(i).(mFields{j}) = nanmean([subD.(mFields{j})]);
                    case 'E'
                        d(i).(mFields{j}) = nanmean([subE.(mFields{j})]);
                end
            end
            if ~ismember(mFields{j},flexors)  % sort muscles by extensors/flexors
                countE = countE+1;
                d(i).extensors(countE) = d(i).(mFields{j});
            else
                countF=countF+1;
                d(i).flexors(countF) = d(i).(mFields{j});
            end
            if ismember(mFields{j},Gs)  %also pull out just the gastrocs
                countGs = countGs+1;
                d(i).Gs(countGs) = d(i).(mFields{j});
            end
        end
    end
    if ~isempty (d(i).rLegMass) % only do this for birds that we have muscle masses for
        d(i).sumExtensors = nansum(d(i).extensors);
        d(i).sumFlexors = nansum(d(i).flexors);
        d(i).sumGs = nansum(d(i).Gs);
    end

end

d = rmfield(d,mFields(6:end));
formatOut = 'dd mmm yyyy';
save('Data\muscleMasses92019.mat','d')

%% import tendon csa csv file

%  csa = table2struct(cssTemp);  Loaded data from Data/CrossSectionalAreasAll43019Simple.xlsx
%     save('TendonCrossAllVals','csa')
% %Add in tendon Data
% load('TendonCSAs2019.mat')

load('Data\TendonCrossAllVals.mat');

cross = struct('bird',{},'tendonCross',{});
c=1;
names = fieldnames(csa);

for i=1:length(names)
    cross(c).bird = str2num(names{i}(2:end));
    cross(c).tendonCross = [csa.(names{i})];
    for j=1:length(d)
        if cross(c).bird == d(j).aniNum           
            d(j).tendonCrossMin = nanmin(cross(c).tendonCross);
            d(j).tendonCrossMax = nanmax(cross(c).tendonCross);
            d(j).tendonCrossAv = nanmean(cross(c).tendonCross);
        end
    end
    c=c+1;
end

%% Add in dxa scan info

opts = spreadsheetImportOptions("NumVariables", 8);

% Specify sheet and range
opts.Sheet = "Sheet1";
opts.DataRange = "A2:H26";

% Specify column names and types
opts.VariableNames = ["aniNum", "AreaCm2", "BmcG", "BmdGcm2", "FatMass", "LeanBMC", "TotalMass", "PerFat"];
opts.VariableTypes = ["double", "double", "double", "double", "double", "double", "double", "double"];

% Import the data
DxaData = readtable("C:\Users\smc288\Box\Suzanne_Cox\SuzanneGF\Botox and Growth\SixMonthFunctionalPaper\Analyses\Data\DxaScanData2017Birds.xlsx", opts, "UseExcel", false);

% Clear temporary variables
clear opts

dxa = table2struct(DxaData);
save('Data\dxaData2017Birds.mat','dxa');

load('Data\dxaData2017Birds.mat')
dxa([dxa.aniNum] == 756).aniNum = 765;
fieldsdxa = fieldnames(dxa);
for i=1:length(d)
    temp = dxa([dxa.aniNum]==d(i).aniNum);
    if ~isempty(temp)
       for j=4:length(fieldsdxa)
           d(i).(fieldsdxa{j}) = temp.(fieldsdxa{j});
       end
    end
end


%% Add in sarcomere lengths
load('C:\Users\smc288\Box\Suzanne_Cox\SuzanneGF\Botox and Growth\SixMonthMorphPaper\MatlabCode\Data\FiberSarcomereLengthsPenn_afterRedo.mat')
fields = fieldnames(fl);
for i=1:length(d)
    temp = fl([fl.bird]==d(i).aniNum);
    if ~isempty(temp)
       for j=4:length(fields)
           d(i).(fields{j}) = temp.(fields{j});
       end
    end
end


%% save data
formatOut = 'dd mmm yyyy';

b792 = d([d.aniNum] == 792);
b800 = d([d.aniNum] == 800);

fields = fieldnames(b800);
for i=1:length(fields)
    if isempty(b792.(fields{i}))
        b792.(fields{i}) =  b800.(fields{i});
    end
end

d([d.aniNum]==792) = b792;

d([d.aniNum]==800) = [];
save(['Data\MuscleTendonFunctionalData',datestr(now, formatOut),'.mat'],'d')

%save

 s2 = rmfield(d,{'rVo2s','sVo2s','extensors','Gs','flexors','FatMass'});
 save(['Data\Vo2JumpData3',datestr(now, formatOut),'.mat'],'s2')
%save(['VO2JumpData',datestr(now, formatOut),'.mat'],'s2')
t2 = struct2table(s2);
writetable(t2,['Data\VO2JumpingMuscleTendonData4R',datestr(now, formatOut),'.csv']);

