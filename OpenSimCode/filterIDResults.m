function [a] = filterIDResults(stoFile,cutoffFreq)
%load GRF mot file
%stoFile = 'B768_S1_inverse_dynamics.sto'
% cutoffFreq = 10;
ID  =importdata(stoFile);
sampleRate = 1/ID.data(2,1);


%apply same filter to both sets of data
[b,a] =  butter(3,cutoffFreq/(300/2));
IDf = ID;
for i=1:size(ID.data,2)
    y=ID.data(:,i);
    if ~sum(~isnan(filtfilt(b,a,y)))
        nanx = isnan(y);
        t = 1:numel(y);
        y(nanx) = interp1(t(~nanx), y(~nanx), t(nanx));
    end
    if (~isempty(y))
        IDf.data(:,i)= filtfilt(b,a,y);
    end
end

%plot a few example columns from each

plot(IDf.data(:,20))
hold on
plot(ID.data(:,20),'s')

plot(IDf.data(:,8))
hold on
plot(ID.data(:,8),'o')


%rewrite GRF Mot and IK mot files
%cutoffFreq of 10 looks good

%Write new IDsto mot file
FullFileName = [stoFile(1:strfind(stoFile,'.')-1),'_f',num2str(cutoffFreq),'hz.mot'];
fid = fopen(FullFileName,'w');
headerLen = size(ID.textdata,1);
for i=1:headerLen-1
     fprintf(fid,'%s\n',ID.textdata{i});
end
%print marker names
%fprintf(fid,'\n\n');
markers = [ID.textdata(headerLen,:)];
headerFormat =[repmat('%s\t',1,length(markers)),'\n'];
fprintf(fid,headerFormat,markers{:});

% Print data
ncol = size(IDf.data,2);
nrows = size(IDf.data,1);
for i= 1:nrows
    for j=1:ncol
        fprintf(fid, '%20.8f\t', IDf.data(i,j));
    end
    fprintf(fid,'\n');
end
fclose(fid);

    

end
