function [ds,numPlates,sampleRate] = forcesFile2mat(fullFileName)
delimiterIn = '\t';
headerlinesIn = 5;
d=importdata(fullFileName,delimiterIn,headerlinesIn);
numplatesTxt = cell2mat(d.textdata(2,1));
numPlates =str2num(numplatesTxt(strfind(numplatesTxt,'=')+1));
samplerateTxt = cell2mat(d.textdata(3,1));
sampleRate = str2num(samplerateTxt(strfind(samplerateTxt,'=')+1:end));
ds = struct('t',{});
ds(1).samples = d.data(:,1);
chan = {'FX','FY','FZ','X','Y','Z','MZ'};
for n=1:numPlates
    for i=1:length(chan)
        y = d.colheaders;
        x1 = [chan{i},num2str(n)];
        substrfind = @(x,y) ~cellfun(@isempty,strfind(y,x));
        logicalArray1 = substrfind(x1,y);
        ds(n).(chan{i})= d.data(:,logicalArray1);
        ds(n).samples = d.data(:,1);
        ds(n).t = (0:sampleRate:(length(d.data(:,1))-1)*sampleRate);
        ds(n).t = ds(n).t/1000;
    end
end

end