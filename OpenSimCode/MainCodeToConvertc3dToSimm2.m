%% First read in all the c3d files into a matlab structure
cd('C:\Users\smc288\Box Sync\SuzanneGF\GuinaModelPaper\OpenSimmStuff\c3d2OSim')
clear
cdBase = pwd;
if ispc
    slash = '\';
else
    slash = '/';
end

d = struct('TrialName',{},'AnalogFrameRate',{},'fileName',{},'filePath',{},'VideoFrameRate',{},'MarkerStruct',{},'AnalogStruct',{});
% 'pt1y',{},'pt2y',{},'pt3y',{},'pt4y',{},'pt5y',{},'pt6y',{},...
% 'pt1z',{},'pt2z',{},'pt3z',{},'pt4z',{},'pt5z',{},'pt6z',{},'toad',{},'hop',{}, 'triggerFrameNum',{});

expfolder = uigetdir ('Pick the folder with the data');
%choose 'B768_c3dFiles' in the folder 'HypoTrackRunning6-18'

tempdirFiles = dir(expfolder);
tempdirFiles(ismember({tempdirFiles.name},{'.','..'})) = [];
tempdirFiles(ismember({tempdirFiles.name},'.DS_Store')) = [];
%topdirFiles = tempdirFiles([tempdirFiles.isdir]);
topdirFiles = tempdirFiles;
idx = 1;
for i = 1:length(topdirFiles)
fileName = tempdirFiles(i).name;
filePath = expfolder;
%[fileName, filePath] = uigetfile('*.c3d', 'Select File', cd, 'MultiSelect', 'on');
fullFileName = strcat(filePath,slash, fileName);
TrialName = regexprep(fileName,{'\.','txt'},{'','txt'});
TrialName = TrialName(1:strfind(TrialName,'c3d')-1);

[MarkerStruct.(TrialName), VideoFrameRate.(TrialName), AnalogStruct.(TrialName),AnalogFrameRate.(TrialName)] ...
    = readc3duwa2Struct_MAC(fullFileName, 'default','y');

end

save('B768TrackRunningDataJuneB2018.mat','MarkerStruct','VideoFrameRate','AnalogStruct','AnalogFrameRate')
    
%% Find ankle joint center and add to struct
clear
load ('B768TrackRunningDataJuneB2018.mat')
% put everything in terms of tibia

% read in knee lateral boundary data
%pull out lat knee boundry data
files = fieldnames(MarkerStruct);
for i=1:length(files)
    if ~isempty(strfind(files{i},'AnkleLat'))
        d = MarkerStruct.([files{i}]);
    end
end
% in the tibia cluster coord system - the translation of the lateral knee
% boundary from the origin
AnkLatBndry =coord_rot_trans(d.TibP, d.TibM, d.TibD,d.V_PntTip);  % here V_PntTip is the virtural marker at the pointer tip
% read in Ankle medial boundary data  
for i=1:length(files)
    if ~isempty(strfind(files{i},'AnkleMed'))
        d = MarkerStruct.([files{i}]);
    end
end
AnkMedBndry =coord_rot_trans(d.TibP, d.TibM, d.TibD,d.V_PntTip);

% plot3(AnkLatBndry(:,1),AnkLatBndry(:,2),AnkLatBndry(:,3))
% hold on
% plot3(AnkMedBndry(:,1),AnkMedBndry(:,2),AnkMedBndry(:,3))

% read in Ankle axis data
for i=1:length(files)
    if ~isempty(strfind(files{i},'AnkleJoint'))
        a = MarkerStruct.([files{i}]);
        disp(i)
    end
end


% define markers (convert ot double percision)
prox1 = double(a.TibP);
mid1 = double(a.TibM);
dist1 = double(a.TibD);
prox2= double(a.TmtP);
mid2 = double(a.TmtM);
dist2 = double(a.TMP);

ankleAxisline = helical2clusters(prox1,mid1,dist1,prox2,mid2,dist2,AnkMedBndry,AnkLatBndry);
[AnkAxisCenter,AnkAxisLine2, AnkerrorMean] = JointCenter(ankleAxisline,AnkLatBndry,AnkMedBndry);

    
% find knee joint center

clearvars -except AnkAxisCenter AnkAxisLine2  AnkerrorMean AnkMedBndry AnkLatBndry
load ('B768TrackRunningDataJuneB2018.mat')
% put everything in terms of tibia

% read in knee lateral boundary data
%pull out lat knee boundry data
files = fieldnames(MarkerStruct);
for i=1:length(files)
    if ~isempty(strfind(files{i},'KneeLat'))
        d = MarkerStruct.([files{i}]);
    end
end
% in the tibia cluster coord system - the translation of the lateral knee
% boundary from the origin
KneeLatBndry =coord_rot_trans(d.TibP, d.TibM, d.TibD,d.V_PntTip);
% read in knee medial boundary data
for i=1:length(files)
    if ~isempty(strfind(files{i},'KneeMed'))
        d = MarkerStruct.([files{i}]);
        disp(i)
    end
end
KneeMedBndry =coord_rot_trans(d.TibP, d.TibM, d.TibD,d.V_PntTip);

% plot3(KneeLatBndry(:,1),KneeLatBndry(:,2),KneeLatBndry(:,3))
% hold on
% plot3(KneeMedBndry(:,1),KneeMedBndry(:,2),KneeMedBndry(:,3))
% read in knee axis data

for i=1:length(files)
    if ~isempty(strfind(files{i},'KneeJoint'))
        a = MarkerStruct.([files{i}]);
        disp(i)
    end
end


% define markers (convert ot double percision)
prox1 = double(a.TibP);
mid1 = double(a.TibM);
dist1 = double(a.TibD);
prox2= double(a.Hip);
mid2 = double(a.FemP);
dist2 = double(a.FemD);

KneeAxisline = helical2clusters(prox1,mid1,dist1,prox2,mid2,dist2,KneeMedBndry,KneeLatBndry);
[KneeAxisCenter,KneeAxisLine2, KneeerrorMean] = JointCenter(KneeAxisline,KneeLatBndry,KneeMedBndry);
% 
% % both knee and ankle joint centers are in local tibia frame of reference
% % define as below
% LTH1 = double(proxTib);
% LTH2 = double(midTib);
% LTH3 = double(distTib);
% originThigh = (LTH1+LTH2+LTH3)/3; % The origin of the local coordinate system.  Usually the distal joint center of the segment, but can be the proximal joint.
% y_axis = (LTH1-LTH3);  % The first defining axis.  Usually the long axis. #smc - use joint edge markers for this
% z_axis = LTH2-originThigh;



% convert tibia points in the walking trial to local tibia coord system - get the transformation matrix
for i=1:length(files)
    if ~isempty(strfind(files{i},'walking'))
        d = MarkerStruct.([files{i}]);
        filename = files{i};
        fileidx = i;
        disp(i)
    end
end
prox1 = double(d.TibP);
mid1 = double(d.TibM);
dist1 = double(d.TibD);
[TibP_L,TibM_L,TibD_L] = coord_rot_trans_3cluster(prox1,mid1,dist1,prox1, mid1, dist1);

% creat an array of knee and ankle centers as long as other pointer data
KneeAxisCenterL = ones(length(d.TibP),3).*KneeAxisCenter';
AnkAxisCenterL = ones(length(d.TibP),3).*AnkAxisCenter';
KneeMed =  ones(length(d.TibP),3).*mean(KneeMedBndry);
KneeLat = ones(length(d.TibP),3).*mean(KneeLatBndry);
AnkMed = ones(length(d.TibP),3).*mean(AnkMedBndry);
AnkLat=  ones(length(d.TibP),3).*mean(AnkLatBndry);

% convert all back into global coord system
d.kneCtr=coord_tr_L_to_Gzc(prox1,mid1,dist1,KneeAxisCenterL);
d.ankCtr=coord_tr_L_to_Gzc(prox1,mid1,dist1,AnkAxisCenterL);

d.kneLat=coord_tr_L_to_Gzc(prox1,mid1,dist1,KneeLat);
d.kneMed=coord_tr_L_to_Gzc(prox1,mid1,dist1,KneeMed);

d.ankLat=coord_tr_L_to_Gzc(prox1,mid1,dist1,AnkLat);
d.ankMed=coord_tr_L_to_Gzc(prox1,mid1,dist1,AnkMed);

save([filename,'MarkersJntCntrsOrigioalRF.mat'],'d')

save('DatawJointCenters.mat','d','filename','VideoFrameRate', 'AnalogFrameRate')

%% Write files
clear
load('DatawJointCenters.mat')
% transform the data to "Model Coordinate system - front of animal is X, right is Z and up is Y
% clearvars -except filename VideoFrameRate AnalogFrameRate
% load ([filename,'MarkersJntCntrsOrigioalRF.mat']);
% data was collected with z up - X in direction animal was running and Y to
% the right
% So just transpose Y and Z
t = [1 ,0,0;0,0,-1;0,1,0];
do = d;  % origional - in lab coordinate frame
fnames = fieldnames(d);
for i=1:length(fnames)
d.(fnames{i}) =  d.(fnames{i})*t;
end

save([filename,'MarkersJntCntrsModelRF.mat'],'d','filename','VideoFrameRate','AnalogFrameRate')

% % Plot 2 check
% markers = fieldnames(d);
% for i=1:length(d.PCaud)
%     for m = 1:length(markers)
%         col = [0.5 m/18 1-m/18];
%         plot3(d.(markers{m})(i,1),d.(markers{m})(i,3),d.(markers{m})(i,2),'*', 'Color',col)
%         hold on
%         xlim([-300,800])
%         ylim([-400,400])
%     end
%     view(0,0)
%     axis equal
%       title(['Frame ',num2str(i)])
%     pause(0.01)
%     hold off
% end



% 
% %calculate angles between joints to set osim model to the same pose
% %draw line between points to visualize angles
% for m = 1:length(markers)
%     col = [0.5 m/14 1-m/14];
%     plot3(dmid.(markers{m})(1),dmid.(markers{m})(2),-dmid.(markers{m})(3),'*', 'Color',col)
%     hold on
%   %  xlim([250,450])
%  %   ylim([-10,400])
%   % zlim([250 450])
% end
% view(0,90)
% xlabel('X')
% ylabel('Y')
% zlabel('Z')
% axis equal
% plot3([dmid.PCaud(1),dmid.Pmid(1)],[dmid.PCaud(2),dmid.Pmid(2)],[dmid.PCaud(3),dmid.Pmid(3)])
% plot3([dmid.Hip(1),dmid.kneCtr(1)],[dmid.Hip(2),dmid.kneCtr(2)],[dmid.Hip(3),dmid.kneCtr(3)])
% plot3([dmid.ankCtr(1),dmid.kneCtr(1)],[dmid.ankCtr(2),dmid.kneCtr(2)],[dmid.ankCtr(3),dmid.kneCtr(3)])
% plot3([dmid.ankCtr(1),dmid.TMP(1)],[dmid.ankCtr(2),dmid.TMP(2)],[dmid.ankCtr(3),dmid.TMP(3)])
% plot3([dmid.ToeP(1),dmid.TMP(1)],[dmid.ToeP(2),dmid.TMP(2)],[dmid.ToeP(3),dmid.TMP(3)])
% plot3([dmid.ToeD(1),dmid.ToeP(1)],[dmid.ToeD(2),dmid.ToeP(2)],[dmid.ToeD(3),dmid.ToeP(3)])

%looks fine - but simm REALLY doesn't like reading in this data
% maybe limitations on rotation of pelvis?
% I had to make it another body and add translation and rotation - I might
% have messed that up

% these are all in mm units - these get converted in writetrc_zmc.m
% Write the TRC file

markers = fieldnames(d);
VideoFrameRatei = VideoFrameRate.(filename);
AnalogFrameRatei = AnalogFrameRate.(filename);
MLabels = fieldnames(d);

%Filter the data
%d2=d;
% datanames = fieldnames(d);
%    [b,a] =  butter(3,60/300);
%    d2=d;
% for j=1:length(datanames)
%     temp =d.(datanames{j});
%     numCords = 3;
%     for m = 1:numCords
%         y= temp(:,m);
%         if ~sum(~isnan(filtfilt(b,a,y)))
%             nanx = isnan(y);
%             t = 1:numel(y);
%             y(nanx) = interp1(t(~nanx), y(~nanx), t(nanx));
%         end
%         if (~isempty(y))
%             d2.(datanames{j})(:,m)= filtfilt(b,a,y);
%         end
%     end
% end


% Write trc files

dT = struct2table(d);
dA = table2array(dT);


save([filename,'wJointMarkers.mat'],'d','VideoFrameRatei','MLabels','dT','dA');

FullFileName = [filename,'ALLJntCntrs.trc'];
writetrc_zmc(dA,MLabels,VideoFrameRatei,FullFileName)

% Write Static Pose TRC file

%mid stance = frame 161
%subset data to just one frame for static pose
dmid = d;
markers = fieldnames(dmid);
for i=1:length(fieldnames(dmid))
    dmid.(markers{i}) = dmid.(markers{i})(161,:);
end

markers = fieldnames(dmid);
VideoFrameRatei = VideoFrameRate.(filename);
MLabels = fieldnames(dmid);
dT = struct2table(dmid);
dA = table2array(dT);
save([filename,'wJointMarkers.mat'],'d','VideoFrameRatei','AnalogFrameRatei','MLabels','dT','dA');

FullFileName = [filename,'Staticf_JntCntrs.trc'];
writetrc_zmc(dA,MLabels,VideoFrameRatei,FullFileName)

% 
% % Break Kinematics into strides and write seporate TRC files
% 
% % break it up into strides
% names = fieldnames(d2);
% s1d = struct('PCaud',{0});
% s2d = struct('PCaud',{0});
% for i = 1:length(names)
%     s1d.(names{i}) = d2.(names{i})(1:217,:);
%     s2d.(names{i}) = d2.(names{i})(189:410,:);
% end
% 
% %checking values
% 
% plot(s2d.TMP(:,1),s2d.TMP(:,2))
% hold on
% plot(s2d.ToeP(:,1),s2d.ToeP(:,2))
% plot(s2d.ToeD(:,1),s2d.ToeD(:,2))
% 
% markers = fieldnames(s1d);
% VideoFrameRatei = VideoFrameRate.(filename);
% MLabels = fieldnames(s1d);
% dT = struct2table(s1d);
% dA = table2array(dT);
% save([filename,'_S1f_wJointMarkers.mat'],'d','VideoFrameRatei','MLabels','dT','dA');
% 
% FullFileName = [filename,'S1f_JntCntrs.trc'];
% writetrc_zmc(dA,MLabels,VideoFrameRatei,FullFileName)
% 
% markers = fieldnames(s2d);
% VideoFrameRatei = VideoFrameRate.(filename);
% MLabels = fieldnames(s2d);
% dT = struct2table(s2d);
% dA = table2array(dT);
% save([filename,'_S2f_wJointMarkers.mat'],'d','VideoFrameRatei','MLabels','dT','dA');
% 
% FullFileName = [filename,'S2f_JntCntrs.trc'];
% writetrc_zmc(dA,MLabels,VideoFrameRatei,FullFileName)


%% Write MOT Files
[a,f,numPlates] = collectForceData;
% for i=1:size(a,2)
%      a2(:,i) = resample(a(:,i),1,3000/300);  %resample the force data at the same framerate as video
% end

%end stride 1 = 0.5640
%start stride 2 0.5640
%end stride 2: 1.2287
% [~,idxend1] = min(abs(0.6480-a(:,1)));
% s1f = a(1:idxend1,:);
% [~,idxstart2] =  min(abs(0.5640-a(:,1)));
% [val,idxend2] =  min(abs(1.2287-a(:,1)));
% AnSampleRate=3000;
% s2f = a(idxstart2:idxend2,:);
% s2f(:,1) = 0:1/AnSampleRate:(length(s2f)-1)/AnSampleRate;
% s1f(:,1) =  0:1/AnSampleRate:(length(s1f)-1)/AnSampleRate;
% writeMot_zmc(s1f(:,2:end),s1f(:,1),[filename,'_S1f_grf.mot'])
% writeMot_zmc(s2f(:,2:end),s2f(:,1),[filename,'_S2f_grf.mot'])
writeMot_zmc(a(:,2:end),a(:,1),[filename,'_All_grf.mot'])

%Check kinematics vs force timing

% plot(sf1.ToeD(:,2))
% hold on
% plot(f(1).FY(1:10:end)*10,'Color','r')
% %plot(f(2).FY(1:10:end)*10)
% plot(f(3).FY(1:10:end)*10,'Color','k')
%plot(f(4).FY(1:10:end)*10)
% externalLoadNames = {'FP_1','FP_2','FP_3','FP_4'};
% AppliedtoBodies = {'phal2_III','phal2_III','phal2_III','phal2_III'};
% GRFFile = [filename,'_grf.mot'];
% MOTFile = 'B768walking10Best2StrikeScaled.mot';
% still working here
%XML = grf2xml_zmc('ExternalLoadNames',externalLoadNames,'AppliedtoBodies',AppliedtoBodies,'GRFFile',GRFFile, 'MOTFile',MOTFile);

%%
% 
% for i=1:length(d.Hip)
%     for m = 1:length(markers)
%         if strcmp(markers{m},'ankCtr')| strcmp(markers{m},'kneCtr')
%             Col = 'r';
%         else
%             Col = 'k';
%         end
%     plot3(d.([markers{m}])(i,1),d.([markers{m}])(i,2),d.([markers{m}])(i,3),'o', 'Color',Col)
% 
%     hold on
%     grid
%     end
%      axis square
%      grid
%         view(0,0)
%     pause(0.1)
%    
% 
%     hold off
% end



%% other code

%This might also help with a bit more description.  This was done for human knee


% given that we have two axes defined, and the centre of the knee joint, we can define a plane perpendicular
% to the Epicondyle axis that intersects the centre of the knee joint. The point that this plane intersects onto 
% the helical axis is defined as the new joint centre, and the helical axis can be redefined in body builder using 
% two new virtual markers.

% The line of the optimal helical axis defined as:
%             x = xo + at
%             y = yo + bt               Eq(1)
%             z = zo + ct
% where xo,yo,zo describe a point along the line (such as Sopt) and <a,b,c> describe the unit vector of the line (helical_axis_vector)
%
%The plane perpendicular to the epicondylar axis and running through the centre of the knee is defined by:
%           a(x-xo)+b(y-yo)+c(z-zo)=0           Eq(2)
% where xo,yo,zo define a point on the plane (knee_centre) and <a,b,c> define the unit vector of the perpendicular plane (epicondylar_axis)
%
% need to join Eq (1) and (2) to solve for t, then substitute this value back into eq (1)

% t = (dot(epicondylar_axis,knee_centre)-dot(epicondylar_axis,Sopt))/dot(epicondylar_axis,helical_axis_vector);
% 
% new_knee_centre_x = Sopt(1)+helical_axis_vector(1)*t;
% new_knee_centre_y = Sopt(2)+helical_axis_vector(2)*t;
% new_knee_centre_z = Sopt(3)+helical_axis_vector(3)*t;
% new_knee_centre=[new_knee_centre_x;new_knee_centre_y;new_knee_centre_z];

% define new virtual markers on the helical axis (make them roughly the width of the knee).
% LHA1=new_knee_centre+helical_axis_vector;
% LHA2=new_knee_centre-helical_axis_vector;
% 
% axisline=[LHA1 LHA2];
% plot3(axisline(1,:),axisline(2,:),axisline(3,:),'b','LineWidth',3)
% scatter3(new_knee_centre(1),new_knee_centre(2),new_knee_centre(3),'bo')
% text(new_knee_centre(1),new_knee_centre(2),new_knee_centre(3),'new knee centre')
% scatter3(LHA1(1),LHA1(2),LHA1(3),'ko')
% text(LHA2(1),LHA2(2),LHA2(3),'LHA 2')
% scatter3(LHA2(1),LHA2(2),LHA2(3),'ko')
% text(LHA1(1),LHA1(2),LHA1(3),'LHA 1')
% 
% title('Left Knee Joint Instantaneous and Optimal Helical Axes')
% LeftKneeCentre=new_knee_centre;
