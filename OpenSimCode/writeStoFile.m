function[] = writeStoFile(ID, FullFileName)
%Write new IDsto mot file
%FullFileName = [stoFile(1:strfind(stoFile,'.')-1),'_f',num2str(cutoffFreq),'hz.mot'];
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
ncol = size(ID.data,2);
nrows = size(ID.data,1);
for i= 1:nrows
    for j=1:ncol
        fprintf(fid, '%20.8f\t', ID.data(i,j));
    end
    fprintf(fid,'\n');
end
fclose(fid);
