function [A] = genMatrix(filename,imgL,imgW,viwN,detN,isflo)
num = str2num(cell2mat(regexp(filename,'\d','match')));
f1 = [filename, '.cof'];
f2 = [filename, '.col'];
f3 = [filename, '.row'];
fid1 = fopen(f1);
fid2 = fopen(f2);
fid3 = fopen(f3);
if strcmp(isflo,'float')
    weg = fread(fid1,num,'float');
else
    weg = fread(fid1,num,'double');
end
colIdx = fread(fid2,num,'int');
rawIdx = fread(fid3,num,'int');
colIdx = colIdx + 1;
rawIdx = rawIdx + 1;
A = sparse(rawIdx,colIdx,weg,viwN*detN,imgL*imgW);
clear colIdx rawIdx weg;
end