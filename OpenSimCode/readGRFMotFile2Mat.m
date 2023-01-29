function [GRF] = readGRFMotFile2Mat(fullFileName)
delimiterIn = '\t';
headerlinesIn = 7;
d=importdata(fullFileName,delimiterIn,headerlinesIn);
GRF =  d.data;
end