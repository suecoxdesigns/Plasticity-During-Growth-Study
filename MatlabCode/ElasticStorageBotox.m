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
    
    modelName = char(model.getName());
    aniNum = str2double(modelfile(regexp(modelfile,'\d')));
    idx = find([d3.AniNum]==aniNum);
    
    ankleRange = [100:5:134];
    tmpRange = [-45:5:-28];
    kneeRange = [-145:5:-130];
    bodyName = 'tibia';
    wrapObjectName = 'ANKLE_4';
    [wrapObject,properties] = getWrapObject(model,bodyName,wrapObjectName);
    angles = d3(idx).maAngles;
    state=model.initSystem;
    jointCoords= model.updCoordinateSet().get('ankle_flexion');
    LG = model.getMuscles().get('LG');
    
    
    %optimize moment arm fit
    
    %initial Guess
    x0 = [wrapObject.get_radius(),osimVec3ToArray(wrapObject.get_translation()),osimVec3ToArray(wrapObject.get_xyz_body_rotation())];
    
    
    radiusMinMax = [x0(1)*0.9,x0(1)*1.3];
    transXMinMax =  sort([x0(2)*0.9,x0(2)*1.3]);
    transYMinMax = sort([x0(3)*0.9,x0(3)*1.3]);
    transZMinMax = [0,0.1];
    
    rotXMinMax = sort([x0(5)*0.9,x0(5)*1.3]);
    rotYMinMax =sort([x0(6)*0.9,x0(6)*1.3]);
    rotZMinMax = [-1,1];
    
    ub = [radiusMinMax(2),transXMinMax(2),transYMinMax(2),transZMinMax(2),rotXMinMax(2),rotYMinMax(2),rotZMinMax(2)];
    lb = [radiusMinMax(1),transXMinMax(1),transYMinMax(1),transZMinMax(1),rotXMinMax(1),rotYMinMax(1),rotZMinMax(1)];
    
    options = optimoptions('fmincon','Display','iter');
    f = @(x0)(calcDiffMa(d3(idx),model,wrapObject,x0));
    try
        [xFopt,fVal] = fmincon(f,x0,[],[],[],[],lb,ub,[],options);
        
        %update wrap object to match optimal values
        wrapObject.set_radius(xFopt(1));
        wrapObject.upd_radius();
        wrapObject.set_xyz_body_rotation(Vec3(xFopt(2),xFopt(3),xFopt(4)));
        wrapObject.upd_xyz_body_rotation();
        wrapObject.set_xyz_body_rotation(Vec3(xFopt(5),xFopt(6),xFopt(7)));
        wrapObject.upd_translation();
        
        %calculate the new moment arm/angle
        count=1;
        for ang = [d3(idx).maAngles]
            model.updCoordinateSet().get('ankle_flexion').setValue(state,ang);
            LgMa(count)=LG.computeMomentArm(state,jointCoords);
            count=count+1;
        end
        diffMa = sum((LgMa+[d3(idx).ma]).^2);
        plot(d3(idx).maAngles,d3(idx).ma)
        hold on
        plot(d3(idx).maAngles,-LgMa,'--')

    catch
        warning('There was a problem. Moving on');
    end
end
save([filename,num2str(i),'.mat'],'data')


    
    
    
    % compute difference between modeled and experimental moment arms
[LgMa,diffMa] = calcDiffMa(d3(idx),model)
    
    
    
    plot(d3(idx).maAngles,d3(idx).ma)
    hold on
    plot(d3(idx).maAngles,-LgMa,'--')
    title([num2str(aniNum),' rms=',num2str(diffMa)])
    
    %calcDiffMa(x0,muscle,jointCoords,state,wrapObject,mas,angleMin,angleMax)
    
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
