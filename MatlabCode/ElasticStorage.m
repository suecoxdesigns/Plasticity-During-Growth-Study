% adjust moment arm to match experimentally measure moment arm across joint
% angles
clear

%load('MatlabCode\Data\MorphData10 Feb 2020.mat');
import org.opensim.modeling.*


simFiles = dir('BotoxBirdModels\*.osim');
load('ElasticDataMay2020.mat')

for i=1:length(simFiles)
      modelfile = simFiles(i).name;
    path = simFiles(i).folder;
    modelFullName = fullfile(path,modelfile);
    shortFileName = modelfile(1:end-5);
    fileOutPath = fullfile(cd,'BotoxBirdModels',[shortFileName,'doneMA.osim']);
    
    model = Model(modelFullName);
    
end


%%


% adjust TSL to match muscle fiber length to measured average length
clear

%load('MatlabCode\Data\MorphData10 Feb 2020.mat');
import org.opensim.modeling.*
simFiles = dir('GuineaModels\GuineaModelsDoneTendonFL\*.osim');
load('MatlabCode\Data\FiberAndSarcomereLengths_afterRedo.mat');
MorphdataFiles = dir('MatlabCode\Data\MorphData*.mat');

[jumk,idx] = max([MorphdataFiles.datenum]);
load(['MatlabCode\Data\',MorphdataFiles(idx).name]);
clear MorphdataFiles idx jumk

load('MatlabCode\Data\MuscleTendonFunctionalData14 Jan 2020.mat')
for i=1:length(simFiles)
    
    modelfile = simFiles(i).name;
    path = simFiles(i).folder;
    modelFullName = fullfile(path,modelfile);
    shortFileName = modelfile(1:end-5);
    fileOutPath = fullfile(cd,'GuineaModels','ModelsAllDataNoMa',[shortFileName,'all.osim']);
    
    model = Model(modelFullName);
    
    state=model.initSystem;
    model.equilibrateMuscles(state);
    modelName = char(model.getName());
    aniNum = str2num(modelName(end-2:end));
    idx = find([d3.aniNum]==aniNum);
    idxd2 =  find([d2.AniNum]==aniNum);
    idxd = find([d.aniNum]==aniNum);
    
    ankleRange = [100:5:134];
    tmpRange = [-45:5:-28];
    kneeRange = [-145:5:-130];
    d3(idx).aniNum = aniNum;
    d3(idx).maxPower = d(idxd).maxPower;
    d3(idx).maxPowerPKg = d(idxd).maxPower/(d(idxd).sumExtensors/1000);
    d3(idx).group = d2(idxd2).group;
    maxPE = 0;
    for a=ankleRange
        for t = tmpRange
            for k=kneeRange
                model.updCoordinateSet().get('hip_flexion').setValue(state,deg2rad(38));
                model.updCoordinateSet().get('knee_flexion').setValue(state,deg2rad(k));
                model.updCoordinateSet().get('ankle_flexion').setValue(state,deg2rad(a));
                model.updCoordinateSet().get('tmp_III_flexion').setValue(state,deg2rad(t));
                
                
                model.updCoordinateSet().get('tmp_III_flexion').getValue(state);
                
                
                
                MG = model.getMuscles().get('MG');
                LG = model.getMuscles().get('LG');
                LG.setActivation(state,100);
                MG.setActivation(state,100);
                state=model.initSystem;
                LG.computeEquilibrium(state);
                MG.computeEquilibrium(state);
               % model.equilibrateMuscles(state);
                
                TendonForce= LG.getTendonForce(state)+ MG.getTendonForce(state);
                [val,idxmin] = min(abs(d2(idxd2).tendonForce-TendonForce));
                strain = d2(idxd2).tendonStrain(idxmin);
                PE = trapz(d2(idxd2).tendonStrain(1:idxmin).*(d2(idxd2).tendonSlackLengthLGmm/1000),d2(idxd2).tendonForce(1:idxmin))
                
                if PE>maxPE
                    d3(idx).PE = PE;
                    d3(idx).PEvals = [k,a,t];
                    maxPE = PE;
                end
            end
        end
    end
     model.disownAllComponents();
end
  save('ElasticData.mat','d3')
formatOut = 'dd mmm yyyy';
Rdata2 = struct2table(d3);
writetable(Rdata2,['ElasticData',datestr(now, formatOut),'.csv']);
%%

clear
%    load('MatlabCode\Data\VO2DataNormalized29 May 2020.mat')
%     load('ElasticData.mat')
%     load('MatlabCode\Data\functionalData.mat')
%     
    
       load('Data\VO2DataNormalized29 May 2020.mat')
    load('Data\ElasticData.mat')
    load('Data\functionalData.mat')
    d3=d2;
    for i=1:length(d3)
        idx = find([df.AniNum]==d3(i).AniNum);
        if isempty(idx)
            continue
        else
            d3(i).maxV = df(idx).maxV;
            d3(i).maxPower = df(idx).maxPower;
            d3(i).work = df(idx).work;
            d3(i).maxV = df(idx).maxV;
            d3(i).peakForce = df(idx).peakForce;
            d3(i).navRvO2b = df(idx).navRvO2b;
            d3(i).mass = df(idx).mass;
            d3(i).treat = df(idx).Treat;
        end
    end

     d3d = d3([d3.treat]=='D')
     d3e = d3([d3.treat] == 'E')
     
    plot([d3e.PE],[d3e.maxV],'o','Color',[0.5 0.5 0.5],'MarkerFaceColor',[0.5 0.5 0.5])
    hold on
plot([d3d.PE],[d3d.maxV],'o','Color','b','MarkerFaceColor','b')
ylabel('Maximum Velocity m/s')
xlabel('Elastic energy storage J')


save('Data\ElasticDataMay2020.mat','d3')
formatOut = 'dd mmm yyyy';
Rdata2 = struct2table(d3);
writetable(Rdata2,['Data\ElasticData',datestr(now, formatOut),'.csv']);
