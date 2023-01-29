% This code consolodates the fiber length information for each bird, muscle
% and muscle location

[path,folder] = fileparts(cd);
dataFolder = fullfile(path,'Data');
%load in calculated fiber length information
%oranize by bird
clear
load('Data\FiberLengthsLgCalculated2.mat')
lg = fl;
clear fl
load('Data\FiberLengthsMgCalculated2.mat')
mg = fl;
clear fl;

%administrative declaration
% 799 MG data fiber lengths is actually 786
% 799 LG data fiber lengths is actually 760
% 790 LG fiber lengths = 786

[lg([lg.bird]==790).bird]= deal(786);
[mg([mg.bird]==799).bird] =deal(786);
[lg([lg.bird]==799).bird] =deal(760);


for i=1:length(mg)
    for j=1:3
        len = 0;
        for k =2:length(mg(i).x(j))
            len = len+ dist3d([mg(i).x(j)])  
            fl = struct('bird',{},'LG',{},'MG',{});
            LG= struct('muscle',{},'loc',{},'avLength',{});
            MG = LG;
            birdNums = unique([lg.bird, mg.bird]);
            MGlocs = {'AS','PS','AD','PD'};
            for i=1:length(birdNums)
            tempLG = lg([lg.bird]==birdNums(i));
            tempMG = mg([mg.bird]==birdNums(i));
            fl(i).bird = birdNums(i);
            for j=1:length(tempLG)
                fl(i).LG(j).muscle = 'LG';
                fl(i).LG(j).loc = tempLG(j).location;
                fl(i).LG(j).avLen = tempLG(j).avLen;
                fl(i).LG(j).sdLen = std(tempLG(j).length);
                fl(i).avLgFl = nanmean([fl(i).LG.avLen]);
                fl(i).sdLgFl = nanstd([fl(i).LG.avLen]);
                switch  fl(i).LG(j).loc
                    case {'AS','PS','AD','PD'}
                        fl(i).LG(j).locN  = find(contains(MGlocs,fl(i).LG(j).loc));
                    case 'P'
                        fl(i).LG(j).locN  =1;
                    case 'M'
                        fl(i).LG(j).locN  =2;
                    case 'D'
                        fl(i).LG(j).locN  =3;
                end
            end
            for j=1:length(tempMG)
                fl(i).MG(j).muscle = 'MG';
                fl(i).MG(j).loc = tempMG(j).location;
                fl(i).MG(j).avLen = tempMG(j).avLen;
                fl(i).MG(j).sdLen = std(tempMG(j).length);
                fl(i).avMgFl = nanmean([fl(i).MG.avLen]);
                fl(i).sdMgFl = nanstd([fl(i).MG.avLen]);
                switch  fl(i).MG(j).loc
                    case {'AS','PS','AD','PD'}
                        fl(i).MG(j).locN  = find(contains(MGlocs,fl(i).MG(j).loc));
                    case 'P'
                        fl(i).MG(j).locN  =1;
                    case 'M'
                        fl(i).MG(j).locN  =2;
                    case 'D'
                        fl(i).MG(j).locN  =3;
                end
            end
            end
        end
    end
end

save('CompliedFLbyBird.mat','fl')

%% read in LG sarcomere lengths into ta mat file

opts = spreadsheetImportOptions("NumVariables", 2);

% Specify sheet and range
opts.Sheet = "Sheet1";
opts.DataRange = "A2:B73";

% Specify column names and types
opts.VariableNames = ["Muscle", "AverageSarcomerelength"];
opts.SelectedVariableNames = ["Muscle", "AverageSarcomerelength"];
opts.VariableTypes = ["string", "double"];
opts = setvaropts(opts, 1, "WhitespaceRule", "preserve");
opts = setvaropts(opts, 1, "EmptyFieldRule", "auto");

% Import the data
LgSarcomereLengthsSmcafterRedo4 = readtable("C:\Users\smc288\Box\Suzanne_Cox\SuzanneGF\Botox and Growth\SixMonthFunctionalPaper\Analyses\Data\LgSarcomereLengthsSmc_afterRedo4.xlsx", opts, "UseExcel", false);

LgSL = table2struct(LgSarcomereLengthsSmcafterRedo4);

for i=1:length(LgSL)
    if LgSL(i).AverageSarcomerelength<0.5
        LgSL(i).AverageSarcomerelength = LgSL(i).AverageSarcomerelength*10^6;
    end
end

save('Data\LgSarcomereLengths_afterRedo.mat','LgSL')

%% Clear temporary variables
clear opts

%% read MG sarcomere lengths into mat file
%% Import data from spreadsheet
% Script for importing data from the following spreadsheet:
%
%    Workbook: C:\Users\smc288\Box\Suzanne_Cox\SuzanneGF\Botox and Growth\SixMonthFunctionalPaper\Analyses\Data\MGSarcomereLengthsSmc_afterRedo.xlsx
%    Worksheet: Sheet1
%
% Auto-generated by MATLAB on 08-Oct-2019 13:42:43

%% Setup the Import Options
opts = spreadsheetImportOptions("NumVariables", 12);

% Specify sheet and range
opts.Sheet = "Sheet1";
opts.DataRange = "A22:L309";

% Specify column names and types
opts.VariableNames = ["Muscle", "PartoftheFiber", "Slideno", "sampleID", "xlengthmm", "hheigthmm", "xh", "btan1xh", "asinbrads", "Lslam", "Sarcomerelengthmm", "AverageSarcomerelength"];
opts.SelectedVariableNames = ["Muscle", "PartoftheFiber", "Slideno", "sampleID", "xlengthmm", "hheigthmm", "xh", "btan1xh", "asinbrads", "Lslam", "Sarcomerelengthmm", "AverageSarcomerelength"];
opts.VariableTypes = ["string", "string", "string", "string", "double", "double", "double", "double", "double", "double", "double", "double"];
opts = setvaropts(opts, [1, 2, 3, 4], "WhitespaceRule", "preserve");
opts = setvaropts(opts, [1, 2, 3, 4], "EmptyFieldRule", "auto");

% Import the data
MGSarcomereLengthsSmcafterRedo = readtable("C:\Users\smc288\Box\Suzanne_Cox\SuzanneGF\Botox and Growth\SixMonthFunctionalPaper\Analyses\Data\MGSarcomereLengthsSmc_afterRedo.xlsx", opts, "UseExcel", false);


% Clear temporary variables
clear opts = readtable("C:\Users\smc288\Box\Suzanne_Cox\SuzanneGF\Botox and Growth\SixMonthFunctionalPaper\Analyses\Data\MGSarcomereLengthsSmc_afterRedo.xlsx", opts, "UseExcel", false);
MgSarLen = table2struct(MGSarcomereLengthsSmcafterRedo);
clear MGSarcomereLengthsSmcafterRedo
count=1;
MgSL = struct('bird',{},'avSL',{});
 MgSL(1).bird = MgSarLen(1).Muscle;
    MgSL(1).avSL = MgSarLen(1).AverageSarcomerelength;
    count=2;
for i=4:3:length(MgSarLen)
    MgSL(count).bird= MgSarLen(i).Muscle;
    MgSL(count).avSL = MgSarLen(i).AverageSarcomerelength;
    count=count+1;
end

for i=1:length(MgSL)
    if MgSL(i).avSL<0.5
        MgSL(i).avSL = MgSL(i).avSL*10^6;
    end
end
save('Data\MgSarcomereLengths_afterRedo.mat','MgSL')
clear temporary variables
clear opts

%% get sarcomere lengths into the fiber lengths file

clear
load('Data\CompliedFLbyBird.mat')
load('Data\MgSarcomereLengths_afterRedo.mat')
load('Data\LgSarcomereLengths_afterRedo.mat')
MGlocs = {'AS','PS','AD','PD'};

%pull out bird number, muscle and location from sarcomere data
for i=1:length(MgSL)
    
    name= MgSL(i).bird{1}(MgSL(i).bird{1}~=' ');
    if isempty(name)
     
       continue
    end
    MgSL(i).name = name;
    MgSL(i).bird = str2num(MgSL(i).name(1:3));
    MgSL(i).musc= MgSL(i).name(4);
    switch MgSL(i).name(7:end)
        case {'AS','PS','AD','PD'}
            MgSL(i).loc =MgSL(i).name(7:end);
            MgSL(i).locN  = find(contains(MGlocs,MgSL(i).loc));
        case 'Prox'
             MgSL(i).loc = 'P';
             MgSL(i).locN =1;
        case 'MID'
            MgSL(i).loc = 'M';
            MgSL(i).locN =2;
        case 'Distal'
             MgSL(i).loc = 'D';
             MgSL(i).locN =3;
    end
end
MgSL(92:96) = [];
for i=1:length(LgSL)
    
     LgSL(i).name= LgSL(i).Muscle{1}(LgSL(i).Muscle{1}~=' ');
    LgSL(i).bird = str2num(LgSL(i).name(1:3));
    if LgSL(i).bird==790
        LgSL(i).bird = 760;
    end
    LgSL(i).musc= LgSL(i).name(4);
    LgSL(i).avSL = LgSL(i).AverageSarcomerelength;
    switch LgSL(i).name(6:end)
        case {'AS','PS','AD','PD'}
            LgSL(i).loc =LgSL(i).name(6:end);
            LgSL(i).locN  = find(contains(MGlocs,LgSL(i).loc));
        case 'Prox'
             LgSL(i).loc = 'P';
             LgSL(i).locN = 1;
        case 'MID'
            LgSL(i).loc = 'M';
            LgSL(i).locN = 2;
        case 'Distal'
             LgSL(i).loc = 'D';
             LgSL(i).locN = 3;
    end
end
LgSL = rmfield(LgSL,{'AverageSarcomerelength','Muscle'});
LgSL = orderfields(LgSL);
MgSL = orderfields(MgSL);
LgTable = struct2table(LgSL);
MgTable = struct2table(MgSL);
%combine LG and MG data
Sl = [MgTable; LgTable];
Sl = table2struct(Sl);
clear i LgSL MgSL


%link with fiber length data
for i=1:length(fl)
   for j=1:length(fl(i).LG)
       temp =sum([Sl.bird]==fl(i).bird & [Sl.musc]=='L' & [Sl.locN] == fl(i).LG(j).locN);
        if temp>0
            fl(i).LG(j).avSL= Sl([Sl.bird]==fl(i).bird & [Sl.musc]=='L' & [Sl.locN] == fl(i).LG(j).locN).avSL;
            fl(i).LG(j).ofl = fl(i).LG(j).avLen*2.36/fl(i).LG(j).avSL;
        end
   end
   for j=1:length(fl(i).MG)
         temp = sum([Sl.bird]==fl(i).bird & [Sl.musc]=='M' & [Sl.locN] == fl(i).MG(j).locN);
         if temp>0
             fl(i).MG(j).avSL = Sl([Sl.bird]==fl(i).bird & [Sl.musc]=='M' & [Sl.locN] == fl(i).MG(j).locN).avSL;
             fl(i).MG(j).ofl = fl(i).MG(j).avLen*2.36/fl(i).MG(j).avSL;
         end
   end
   if (~isempty(fl(i).LG) && isfield(fl(i).LG,'ofl'))
       fl(i).LgAvOfl = nanmean([fl(i).LG.ofl]);
   end
   if (~isempty(fl(i).MG)&& isfield(fl(i).MG,'ofl'))
       fl(i).MgAvOfl = nanmean([fl(i).MG.ofl]);
   end
end

save('Data\FiberAndSarcomereLengths_afterRedo.mat','fl')

% OFL = FL*2.36/SL;
% Optimal fascicle length (L 0 ) was calculated
% multiplying the length of the fascicle by the ratio of optimal sarcomere length of guinea fowl
% muscle (2.36 �m; Carr et al., 2011b) to the average measured sarcomere length. L 0 were also
% normalized to total limb length to assess whether adaptations occur independent of a general size
% effect, since longer muscles can simply be associated with longer bones (Legerlotz et al., 2016)

% average over the different regions - 

%% read in mtu length data
opts = spreadsheetImportOptions("NumVariables", 6);

% Specify sheet and range
opts.Sheet = "Sheet1";
opts.DataRange = "A2:F526";

% Specify column names and types
opts.VariableNames = ["birdNum", "musc", "X", "Y", "Z", "Point"];
opts.VariableTypes = ["double", "categorical", "double", "double", "double", "categorical"];
opts = setvaropts(opts, [2, 6], "EmptyFieldRule", "auto");

% Import the data
tbl = readtable("C:\Users\smc288\Box\Suzanne_Cox\SuzanneGF\Botox and Growth\SixMonthFunctionalPaper\Analyses\Data\MsucleTendonPathDataCompiled.xlsx", opts, "UseExcel", false);

% Convert to output type
birdNum = tbl.birdNum;
musc = tbl.musc;
X = tbl.X;
Y = tbl.Y;
Z = tbl.Z;
Point = tbl.Point;

% Clear temporary variables
clear opts tbl

%% calculate mtu lengths
%get muscle names all in capital
musc(musc=='lg') = deal('LG');
musc(musc=='mg') = deal('MG');


birds = unique(birdNum);
for i=1:length(birds)
    for mu = {'LG','MG'}
        xyz = [X(birdNum==birds(i) & musc==mu),Y(birdNum==birds(i) & musc==mu),Z(birdNum==birds(i) & musc==mu)];
        mtuLen = xyz;
        % Compute distance between each digitised point
        for k = 1:length(mtuLen)-1
            tmp=[mtuLen(k+1,:) - mtuLen(k,:)].^2;
            dist(k,:)=[sum(tmp')'].^0.5;
        end
        
        pointsSub = Point(birdNum==birds(i)& musc==mu);
        musXyz = xyz(1:find(pointsSub == 'End of muscle belly'),:);
       for k = 1:length(musXyz)-1
            tmp=[musXyz(k+1,:) - musXyz(k,:)].^2;
            distMusc(k,:)=[sum(tmp')'].^0.5;
        end
        % Compute total length
        if sum(ismember([fl.bird],birds(i)))
        fl([fl.bird]==birds(i)).(char(strcat('mtuLeng',mu))) = sum(dist);
        fl([fl.bird]==birds(i)).(char(strcat(mu,'Length'))) = sum(distMusc);
        else
            last = length(fl);
            fl(last+1).bird = birds(i);
            fl(last+1).(char(strcat('mtuLeng',mu))) = sum(dist);
            fl(last+1).(char(strcat(mu,'Length'))) = sum(distMusc);
        end
    end
end

fl([fl.bird]==781).bird = 785;
save('Data\FiberAndSarcomereLengths.mat','fl')


%% read in pennation angle data
folder = 'C:\Users\smc288\Box\Suzanne_Cox\SuzanneGF\Botox and Growth\SixMonthFunctionalPaper\Analyses\Data\Pennation\UpdatedPennationData';
files = dir(folder);
filenames = {files([files.isdir]~=1).name};
for i=1:length(filenames)
    birdNum = str2num(filenames{i}(1:3));
    if birdNum==790
        birdNum = 760;
    end
    if birdNum==786
        i;
    end
    if birdNum == 798
        i;
    end
    underscrs = strfind(filenames{i},'_');
    musc = lower(filenames{i}(5:6));
    if strcmpi(musc,'MG')
        musc = [musc,filenames{i}(8)];
    end
    temp =csvread(fullfile(folder,filenames{i}),1,5);
    avPenAngle = nanmean(temp);
    if ismember(birdNum,[fl.bird])
        fl([fl.bird]==birdNum).([musc,'PenAngle']) = avPenAngle;
    else
        fl(end+1).bird = birdNum;
        fl(end).([musc,'PenAngle']) = avPenAngle;
    end
end

save('Data\FiberSarcomereLengthsPenn_afterRedo.mat','fl')


% %% combine with jumping muscle ect data  Vo2JumpData2119.mat
% clear
% load('FiberSarcomereLengthsPenn.mat')
% load('Vo2JumpData2119.mat')
% fields = fieldnames(fl)
% for i=1:length(s2)
%     temp = fl([fl.bird]==s2(i).aniNum);
%     if ~isempty(temp)
%        for j=4:length(fields)
%            s2(i).(fields{j}) = temp.(fields{j});
%        end
%     end
% end
% 
% save('Vo2JumpData22019.mat','s2')
% 
% formatOut = 'dd mmm yyyy';
% %save(['VO2JumpData',datestr(now, formatOut),'.mat'],'s2')
% t2 = struct2table(s2);
% writetable(t2,['VO2JumpingData4R',datestr(now, formatOut),'.csv']);